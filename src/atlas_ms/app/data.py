"""
Reading a processed project for the app.

Everything here is quick to load: Parquet tables, and the MGF read with a
small parser (``atlas_ms.mgf``; matchms would add about 10 s to the app start). Chromatograms
are not stored by the pipeline: they are extracted from the mzML files when a
feature is selected (``preprocessing.msdata.RunReader``).
"""

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from atlas_ms.mgf import read_mgf
from atlas_ms.project import Project


@dataclass
class Results:
    """The results of a project, as the app uses them."""

    features: pd.DataFrame  # features.parquet + network columns (NaN for features without MS2) + best annotation
    quant: pd.DataFrame  # intensities, indexed by feature_id, one column per sample
    gap_filled: pd.DataFrame | None  # same shape as quant, True where gap-filled
    edges: pd.DataFrame  # network edges (empty before the network exists)
    spectra: dict[int, tuple[np.ndarray, np.ndarray]]  # feature_id -> (m/z, intensity) of its MS2
    candidates: pd.DataFrame  # every annotation candidate (annotations/candidates.parquet; empty before)


# Best-annotation columns (annotations/best.parquet) added to the feature
# table, renamed where the name alone would be unclear.
BEST_COLUMNS = {
    "name": "annotation", "level": "level", "label": "annotation_label", "source": "annotation_source",
    "lipid_class": "lipid_class", "lipid_name": "lipid_name", "lipid_level": "lipid_level",
    "flags": "annotation_flags", "family_class": "family_class", "family_class_score": "family_class_score",
}


def load_results(project: Project) -> Results | None:
    """The project's results, or None if the pipeline has not produced them yet."""
    results = project.results_dir
    if not (results / "features.parquet").exists():
        return None
    features = pd.read_parquet(results / "features.parquet")
    nodes_file = results / "network" / "nodes.parquet"
    if nodes_file.exists():
        features = features.merge(pd.read_parquet(nodes_file), on="feature_id", how="left")
    best_file = results / "annotations" / "best.parquet"
    if best_file.exists():
        best = pd.read_parquet(best_file)[["feature_id", *BEST_COLUMNS]].rename(columns=BEST_COLUMNS)
        features = features.merge(best, on="feature_id", how="left")
    else:  # not annotated yet
        features = features.assign(**{column: np.nan for column in BEST_COLUMNS.values()})
    # Text columns: "" rather than missing, so tables and tooltips show blanks.
    text = [column for column in BEST_COLUMNS.values() if column != "family_class_score"]
    features[text] = features[text].fillna("").astype(str)
    candidates_file = results / "annotations" / "candidates.parquet"
    edges_file = results / "network" / "edges.parquet"
    edges = pd.read_parquet(edges_file) if edges_file.exists() else pd.DataFrame(
        columns=["source", "target", "edge_type", "score", "matched_peaks", "mz_delta"])
    gap_filled_file = results / "quant_gap_filled.parquet"
    mgf = results / "gnps" / "ms2_spectra.mgf"
    return Results(
        features=features,
        quant=pd.read_parquet(results / "quant.parquet").set_index("feature_id"),
        gap_filled=pd.read_parquet(gap_filled_file).set_index("feature_id") if gap_filled_file.exists() else None,
        edges=edges,
        spectra=read_mgf(mgf) if mgf.exists() else {},
        candidates=pd.read_parquet(candidates_file) if candidates_file.exists() else pd.DataFrame(),
    )


def log_intensity(values: pd.Series) -> pd.Series:
    """log10 of intensities, with 0 for missing or non-positive values."""
    return values.where(values > 0).apply(lambda v: math.log10(v) if v == v else 0.0)


def find_fragments(
    spectra: dict[int, tuple[np.ndarray, np.ndarray]],
    masses: list[float],
    tolerance: float,
    precursors: dict[int, float] | None = None,
) -> list[int]:
    """
    Feature ids whose MS2 spectrum contains *every* one of ``masses``
    (within ``tolerance``, in Da).

    With ``precursors`` (feature id -> precursor m/z), ``masses`` are neutral
    losses instead: a loss L is present when a fragment sits at
    precursor - L. Class-specific losses are how many lipids are recognised
    in positive mode (e.g. 141.019 for PE, 185.009 for PS).
    """
    found = []
    for feature_id, (mz, _) in spectra.items():
        values = mz if precursors is None else precursors[feature_id] - mz
        if len(values) and all(np.abs(values - mass).min() <= tolerance for mass in masses):
            found.append(feature_id)
    return found
