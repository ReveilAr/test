"""
Harmonization: all annotation sources -> confidence levels, flags, and one
best annotation per feature.

1. Lipid names from every source are read with Goslin (``lipids.parse_lipid``),
   so that "PC(16:0/18:1(9Z))" from a library and "PC 34:1" from the rules
   are recognised as the same class and species.
2. Levels: each candidate keeps its source's proposed Schymanski level,
   except where a rule caps it. Harmonization may lower confidence, never
   raise it:
   * in-silico library spectra, rule-based lipids, CSI:FingerID and the
     other structure predictions: at most level 3;
   * level 1 needs a retention time match: without one, 2a.
3. Best annotation per feature: the most confident level, then the source
   priority (``source_priority``), then the score.
4. Flags (never a reason to drop a candidate):
   * ``blank``: the feature is as intense in blanks as in samples;
   * ``lipid class disagreement``: the sources' top candidates (level 3 or
     better) name different lipid classes;
   * ``ambiguous``: the best source has another, different candidate with
     the same level and score;
   * ``RT off class trend``: in reversed-phase LC, the retention time of a
     lipid class rises with chain length and falls with unsaturation. Each
     class with enough confident annotations gets a linear fit
     RT ~ carbons + double bonds (the equivalent carbon number idea), robust
     to outliers; best annotations too far from their class's line are
     flagged.
5. Family consensus: within each molecular family, the most frequent lipid
   class of the confident best annotations, and the fraction of annotated
   nodes that agree (``family_class``, ``family_class_score``). This is the
   family-level class propagation of MolNetEnhancer.
"""

import logging
import re

import numpy as np
import pandas as pd

from atlas_ms.annotation.lipids import parse_lipid
from atlas_ms.annotation.schema import LEVEL_RANK, empty_table, worst_level
from atlas_ms.config import HarmonizationSettings

log = logging.getLogger(__name__)

# Sources whose candidates are at most level 3 (structure predictions and
# class-level evidence). Library hits are capped by their kind instead.
CAPPED_AT_3 = ("lipid_rules", "sirius:csi", "sirius:canopus", "sirius:elgordo", "ms2query")

BEST_COLUMNS = ["feature_id", "name", "level", "label", "source", "score", "formula", "adduct", "smiles",
                "inchikey", "lipid_class", "lipid_name", "lipid_species", "lipid_level", "n_candidates",
                "flags", "family_class", "family_class_score"]


def final_level(candidate) -> str:
    """The candidate's level after the caps (never more confident than proposed)."""
    level = candidate.proposed_level
    if candidate.library_kind == "in_silico" or candidate.source.startswith(CAPPED_AT_3):
        level = worst_level(level, "3")
    if level == "1" and not np.isfinite(candidate.rt_error_s):
        level = "2a"  # no retention time: not confirmed by a standard
    return level


def add_lipid_names(candidates: pd.DataFrame) -> pd.DataFrame:
    """
    Fill lipid class, normalised name, species and structural level for
    every candidate whose name Goslin reads as a lipid (library hits named
    in lipid shorthand, for example). ``lipid_species`` is the sum
    composition, used to compare sources.
    """
    parsed = [parse_lipid(name) if name else None for name in candidates["name"]]
    candidates = candidates.copy()
    for column, key in (("lipid_class", "lipid_class"), ("lipid_name", "name"), ("lipid_level", "level")):
        from_goslin = [p[key] if p else "" for p in parsed]
        # Keep what the source already said (the lipid rules fill these themselves).
        candidates[column] = candidates[column].where(candidates[column] != "", from_goslin)
    candidates["lipid_species"] = [p["species"] if p else "" for p in parsed]
    return candidates


def source_family(source: str) -> str:
    """"library:inhouse" -> "library" (the name used in source_priority)."""
    return source.split(":", 1)[0]


def choose_best(candidates: pd.DataFrame, settings: HarmonizationSettings) -> pd.DataFrame:
    """One row per feature: its most confident candidate."""
    priority = {name: index for index, name in enumerate(settings.source_priority)}
    ranked = candidates.assign(
        _level=candidates["level"].map(LEVEL_RANK),
        _priority=candidates["source"].map(lambda s: priority.get(source_family(s), len(priority))),
    ).sort_values(["feature_id", "_level", "_priority", "score"], ascending=[True, True, True, False])
    best = ranked.drop_duplicates("feature_id").drop(columns=["_level", "_priority"])
    best["n_candidates"] = best["feature_id"].map(candidates.groupby("feature_id").size())
    return best.set_index("feature_id")


def feature_flags(candidates: pd.DataFrame, best: pd.DataFrame, blank_ids: set) -> dict[int, list[str]]:
    """Flags per feature id (see the module docstring), except the RT trend."""
    flags = {feature_id: [] for feature_id in best.index}
    for feature_id in blank_ids & set(flags):
        flags[feature_id].append("blank")
    confident = candidates[(candidates["level"].map(LEVEL_RANK) <= LEVEL_RANK["3"]) & (candidates["rank"] == 1)]
    classes = confident[confident["lipid_class"] != ""].groupby("feature_id")["lipid_class"].nunique()
    for feature_id in classes.index[classes > 1]:
        flags[feature_id].append("lipid class disagreement")
    for feature_id, group in candidates.groupby("feature_id"):
        top = best.loc[feature_id]
        rivals = group[(group["source"] == top["source"]) & (group["level"] == top["level"])
                       & (group["score"] >= top["score"] - 1e-9) & (group["name"] != top["name"])]
        if not rivals.empty:
            flags[feature_id].append("ambiguous")
    return flags


_COMPOSITION = re.compile(r" (?:O-)?(\d+):(\d+)")


def deleted_residuals(design: np.ndarray, rt: np.ndarray) -> np.ndarray:
    """
    For each point, its error when predicted by a fit *without* it:
    residual / (1 - leverage). A point far from the others (high leverage)
    pulls an ordinary fit towards itself and hides its own error; this
    measure does not.
    """
    coefficients, *_ = np.linalg.lstsq(design, rt, rcond=None)
    leverage = np.einsum("ij,jk,ik->i", design, np.linalg.pinv(design.T @ design), design)
    return (rt - design @ coefficients) / np.clip(1.0 - leverage, 1e-6, None)


def rt_outliers(best: pd.DataFrame, rt_min: pd.Series, settings: HarmonizationSettings) -> set:
    """
    Feature ids whose best lipid annotation is off its class's RT trend.
    Per class, with at least ``rt_model_min_points`` confident species, a
    least-squares line RT (min) = a + b * carbons + c * double bonds. The
    point with the largest deleted residual is set aside while that residual
    exceeds ``rt_outlier_min`` (and enough points remain), then everything
    is judged against the line of the remaining points.
    """
    lipids = best[(best["lipid_species"] != "") & (best["level"].map(LEVEL_RANK) <= LEVEL_RANK["3"])]
    outliers = set()
    for _, group in lipids.groupby("lipid_class"):
        if len(group) < settings.rt_model_min_points:
            continue
        composition = group["lipid_species"].str.extract(_COMPOSITION).astype(float)
        design = np.column_stack([np.ones(len(group)), composition[0], composition[1]])
        rt = rt_min.reindex(group.index).to_numpy()
        keep = np.ones(len(group), dtype=bool)
        while True:
            deleted = np.zeros(len(group))
            deleted[keep] = deleted_residuals(design[keep], rt[keep])
            worst = np.argmax(np.abs(deleted))
            if abs(deleted[worst]) <= settings.rt_outlier_min or keep.sum() <= settings.rt_model_min_points:
                break
            keep[worst] = False
        coefficients, *_ = np.linalg.lstsq(design[keep], rt[keep], rcond=None)
        off = np.abs(rt - design @ coefficients) > settings.rt_outlier_min  # set-aside points vs the clean line
        off[keep] = np.abs(deleted[keep]) > settings.rt_outlier_min  # the others: their own deleted residual
        outliers |= set(group.index[off])
    return outliers


def family_consensus(best: pd.DataFrame, nodes: pd.DataFrame) -> pd.DataFrame:
    """
    Per molecular family: the most frequent lipid class among the confident
    (level 3 or better) best annotations, and the fraction of those
    annotations that agree. Returned per feature (every node of the family).
    """
    families = nodes.loc[nodes["family"] > 0, ["feature_id", "family"]].set_index("feature_id")["family"]
    annotated = best.reindex(families.index)
    confident = annotated[annotated["lipid_class"].fillna("") != ""]
    confident = confident[confident["level"].map(LEVEL_RANK) <= LEVEL_RANK["3"]]
    rows = []
    for family, group in confident.groupby(families.reindex(confident.index)):
        counts = group["lipid_class"].value_counts()
        rows.append({"family": family, "family_class": counts.index[0],
                     "family_class_score": counts.iloc[0] / counts.sum()})
    consensus = pd.DataFrame(rows, columns=["family", "family_class", "family_class_score"])
    return families.rename("family").reset_index().merge(consensus, on="family").drop(columns="family")


def harmonize(
    candidates: pd.DataFrame, features: pd.DataFrame, nodes: pd.DataFrame, settings: HarmonizationSettings,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    All sources' candidates -> (candidates with final levels, best
    annotation per feature). ``best`` has one row per feature of the
    feature table, empty where nothing was found.
    """
    candidates = add_lipid_names(candidates)
    candidates["level"] = [final_level(c) for c in candidates.itertuples()]
    all_ids = pd.Index(features["feature_id"], name="feature_id")
    if candidates.empty:
        best = pd.DataFrame(index=all_ids, columns=[c for c in BEST_COLUMNS if c != "feature_id"])
        return candidates, best.reset_index()

    best = choose_best(candidates, settings)
    blank = set(nodes.loc[nodes["qc"] == "present in blanks", "feature_id"]) if "qc" in nodes else set()
    flags = feature_flags(candidates, best, blank)
    rt_min = features.set_index("feature_id")["rt"] / 60.0
    for feature_id in rt_outliers(best, rt_min, settings):
        flags[feature_id].append("RT off class trend")
    best["flags"] = [";".join(flags[i]) for i in best.index]
    best["label"] = [
        f"L{row.level} · {row.name}" + (f" ({row.lipid_level})" if row.lipid_level else "")
        for row in best.itertuples()
    ]
    best = best.reindex(all_ids).reset_index()
    best = best.merge(family_consensus(best.set_index("feature_id"), nodes), on="feature_id", how="left")
    log.info("Best annotations: %d features; by level: %s", best["level"].notna().sum(),
             best["level"].value_counts().sort_index().to_dict())
    return candidates, best[BEST_COLUMNS]


def run_harmonize(annotation_files: list, features_file, nodes_file, candidates_out, best_out,
                  settings: HarmonizationSettings) -> None:
    """File-level entry point of the ``harmonize`` rule."""
    tables = [pd.read_parquet(path) for path in annotation_files]
    tables = [table for table in tables if not table.empty]
    candidates = pd.concat(tables, ignore_index=True) if tables else empty_table()
    candidates, best = harmonize(candidates, pd.read_parquet(features_file), pd.read_parquet(nodes_file), settings)
    candidates.to_parquet(candidates_out, index=False)
    best.to_parquet(best_out, index=False)
