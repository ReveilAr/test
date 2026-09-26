"""
Rule-based lipid annotation (positive ion mode).

Standard tools (library search, MS2Query, CSI:FingerID) do poorly on lipids,
whose MS2 spectra are sparse and dominated by a few class-specific ions. This
module reads them the way a lipidomics analyst does, with the rules of
``atlas_ms/presets/lipid_rules.yaml``:

1. Species list: for every class and adduct of the rules, every species in
   the class's carbon / double-bond ranges, with the m/z of its ion. The
   formula of a species is derived from the class's reference species (one
   more carbon: + CH2, one more double bond: - H2).
2. Precursor: a feature's m/z is matched to the species list
   (``precursor_tolerance_ppm``). Several classes can match the same m/z
   (e.g. PC 34:1 and PE 37:1 share a formula).
3. Class evidence: the diagnostic fragments and neutral losses of the class
   and adduct must be in the MS2 spectrum (e.g. m/z 184.0733 for PC
   [M+H]+, a loss of 141.0191 for PE [M+H]+). This separates the isobars.
4. Chains: when fragments reveal the chains, the species (``TG 52:2``)
   becomes a molecular species (``TG 16:0_18:1_18:1``): fatty acid losses of
   glycerolipids, sphingoid base ions of sphingolipids (``Cer 18:1;O2/16:0``).
   A composition is only given when it is the single one explaining the
   observed losses.

Score: the fraction of the rule's ions (and chain evidence) found. A matched
rule proposes Schymanski level 3 (tentative candidate). Features where no
rule matched (including those without MS2) can get the species matching
their m/z only, as level 5 suggestions.

The structural level (species, molecular species, sn-position...) and the
normalised name come from Goslin (pygoslin), the reference lipid shorthand
parser, so that names from all sources can be compared.
"""

import logging
import re
from importlib import resources
from itertools import combinations_with_replacement
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from atlas_ms.annotation.schema import candidate_table, to_json
from atlas_ms.config import LipidSettings
from atlas_ms.mgf import read_mgf

log = logging.getLogger(__name__)

# ---- Masses ----------------------------------------------------------------

# Monoisotopic masses of the elements used by lipids and adducts.
ELEMENT_MASS = {
    "C": 12.0, "H": 1.00782503207, "N": 14.0030740048, "O": 15.99491461956,
    "P": 30.97376163, "S": 31.97207100, "Na": 22.9897692809, "K": 38.9637064864,
}
ELECTRON = 0.00054857990946

# What each adduct adds to the neutral molecule (one positive charge).
ADDUCTS = {
    "[M+H]+": {"H": 1},
    "[M+Na]+": {"Na": 1},
    "[M+NH4]+": {"N": 1, "H": 4},
    "[M+K]+": {"K": 1},
}


def parse_formula(text: str) -> dict[str, int]:
    """"C42H82NO8P" -> {"C": 42, "H": 82, "N": 1, "O": 8, "P": 1}."""
    counts = {}
    for element, number in re.findall(r"([A-Z][a-z]?)(\d*)", text):
        counts[element] = counts.get(element, 0) + (int(number) if number else 1)
    return counts


def format_formula(counts: dict[str, int]) -> str:
    """Hill order: C, H, then the other elements alphabetically."""
    order = ["C", "H"] + sorted(e for e in counts if e not in ("C", "H"))
    return "".join(f"{e}{counts[e] if counts[e] != 1 else ''}" for e in order if counts.get(e))


def formula_mass(counts: dict[str, int]) -> float:
    return sum(ELEMENT_MASS[element] * number for element, number in counts.items())


H2O = formula_mass({"H": 2, "O": 1})
NH3 = formula_mass({"N": 1, "H": 3})
PROTON = ELEMENT_MASS["H"] - ELECTRON


def ion_mz(neutral_mass: float, adduct: str) -> float:
    """m/z of the singly charged ion of a neutral mass."""
    return neutral_mass + formula_mass(ADDUCTS[adduct]) - ELECTRON


# Fatty acids that can be lost from a glycerolipid (C8-C28, 0-6 double bonds),
# as neutral acids C(c)H(2c-2d)O2.
FATTY_ACIDS = [(c, d) for c in range(8, 29) for d in range(0, 7) if d <= c // 2 - 1]
FATTY_ACID_MASS = np.array([formula_mass({"C": c, "H": 2 * c - 2 * d, "O": 2}) for c, d in FATTY_ACIDS])

# Sphingoid bases with two hydroxyls (d-bases, e.g. sphingosine 18:1;O2,
# C18H37NO2), C14-C22 with 0-2 double bonds. They show as [base+H-H2O]+ and
# [base+H-2H2O]+ (m/z 282.2791 and 264.2686 for 18:1;O2).
SPHINGOID_BASES = [(c, d) for c in range(14, 23) for d in range(0, 3)]
_BASE_MASS = np.array([formula_mass({"C": c, "H": 2 * c + 3 - 2 * d, "N": 1, "O": 2}) for c, d in SPHINGOID_BASES])
SPHINGOID_ION_MZ = np.stack([_BASE_MASS + PROTON - H2O, _BASE_MASS + PROTON - 2 * H2O], axis=1)

# Diacylglycerol-like ions [DG+H-H2O]+ (DG c:d = C(c+3)H(2c-2d+4)O5), e.g.
# m/z 577.5190 for DG 34:1: the phosphatidyl halves of cardiolipins.
DG_COMPOSITIONS = [(c, d) for c in range(24, 45) for d in range(0, 9)]
DG_ION_MZ = np.array([formula_mass({"C": c + 3, "H": 2 * c - 2 * d + 4, "O": 5}) + PROTON - H2O
                      for c, d in DG_COMPOSITIONS])


# ---- Rules and species -----------------------------------------------------

def rules_path(rules_file: str = "") -> str:
    """The rules file: the user's own, or the built-in one."""
    return rules_file or str(resources.files("atlas_ms") / "presets" / "lipid_rules.yaml")


def load_rules(path: str | Path) -> dict:
    """Lipid classes from a rules file (see the comments of lipid_rules.yaml)."""
    return yaml.safe_load(Path(path).read_text())["classes"]


def species_table(rules: dict) -> pd.DataFrame:
    """
    Every species of every class and adduct, with its ion m/z, sorted by
    m/z. Columns: lipid_class, adduct, carbons, double_bonds, name,
    formula, mz.
    """
    rows = []
    for lipid_class, rule in rules.items():
        reference = parse_formula(rule["formula"])
        ref_c, ref_d = rule["reference"]
        for c in range(rule["carbons"][0], rule["carbons"][1] + 1):
            for d in range(rule["double_bonds"][0], rule["double_bonds"][1] + 1):
                if d > c // 2:  # at most one double bond per two carbons
                    continue
                counts = dict(reference, C=reference["C"] + (c - ref_c),
                              H=reference["H"] + 2 * (c - ref_c) - 2 * (d - ref_d))
                mass = formula_mass(counts)
                for adduct in rule["adducts"]:
                    rows.append({"lipid_class": lipid_class, "adduct": adduct, "carbons": c, "double_bonds": d,
                                 "name": rule["name"].format(c=c, d=d), "formula": format_formula(counts),
                                 "mz": ion_mz(mass, adduct)})
    return pd.DataFrame(rows).sort_values("mz", ignore_index=True)


# ---- Evidence in a spectrum ------------------------------------------------

def found(peaks_mz: np.ndarray, targets, tolerance: float) -> np.ndarray:
    """For each target m/z: is there a peak within ``tolerance``?"""
    targets = np.atleast_1d(np.asarray(targets, dtype=float))
    if len(peaks_mz) == 0:
        return np.zeros(targets.shape, dtype=bool)
    return (np.abs(peaks_mz[:, None] - targets.reshape(1, -1)) <= tolerance).any(axis=0).reshape(targets.shape)


def chain_name(c: int, d: int) -> str:
    return f"{c}:{d}"


def combine_chains(observed: list, n_chains: int, total: tuple[int, int]):
    """
    The single way of writing the total composition as ``n_chains`` chains
    that includes every observed chain (repeats allowed), either from the
    observed chains alone or with one more chain making up the rest (a chain
    whose loss was not seen). None when there is no way, or several.
    """
    def complete(combo):
        return all(chain in combo for chain in observed)

    options = set()
    for combo in combinations_with_replacement(sorted(observed), n_chains):
        if (sum(c for c, _ in combo), sum(d for _, d in combo)) == total and complete(combo):
            options.add(combo)
    if not options:
        for combo in combinations_with_replacement(sorted(observed), n_chains - 1):
            rest = (total[0] - sum(c for c, _ in combo), total[1] - sum(d for _, d in combo))
            if rest[0] >= 2 and 0 <= rest[1] <= rest[0] // 2 and complete(combo + (rest,)):
                options.add(tuple(sorted(combo + (rest,))))
    return list(options.pop()) if len(options) == 1 else None


def evaluate(rule: dict, species: pd.Series, precursor: float, peaks_mz: np.ndarray, tolerance: float) -> dict:
    """
    Check one species' class rule against a spectrum (already reduced to
    the peaks above the relative-intensity threshold). Returns the evidence:
    matched (all required evidence found), found / missing ions, the name
    (with chains when known), and the score.
    """
    found_items, missing, required_ok = [], [], True
    for key, required, loss in (("fragments", True, False), ("optional_fragments", False, False),
                                ("neutral_losses", True, True), ("optional_neutral_losses", False, True)):
        values = rule.get(key, [])
        targets = [precursor - v for v in values] if loss else values
        for value, hit in zip(values, found(peaks_mz, targets, tolerance) if values else []):
            label = f"NL {value:.4f}" if loss else f"{value:.4f}"
            (found_items if hit else missing).append(label)
            required_ok &= bool(hit) or not required

    if rule.get("dg_ions"):
        hits = found(peaks_mz, DG_ION_MZ, tolerance)
        dg = [chain_name(*composition) for composition, hit in zip(DG_COMPOSITIONS, hits) if hit]
        (found_items if dg else missing).append("DG-like ion" + (f" ({', '.join(dg)})" if dg else ""))
        required_ok &= bool(dg)

    name = species["name"]
    method = rule.get("chains")
    if method:
        total = (int(species["carbons"]), int(species["double_bonds"]))
        chains = _chains(method, species, precursor, peaks_mz, tolerance, total)
        if chains:
            name = chains["name"]
            found_items.append(chains["evidence"])
        else:
            missing.append("chain evidence")
        required_ok &= chains is not None or not rule.get("chains_required", False)

    n_items = len(found_items) + len(missing)
    return {"matched": required_ok, "name": name, "found": found_items, "missing": missing,
            "score": len(found_items) / n_items if n_items else 0.0}


def _chains(method: str, species: pd.Series, precursor: float, peaks_mz: np.ndarray, tolerance: float, total):
    """Chain composition from the fragments, as {"name": ..., "evidence": ...}, or None."""
    lipid_class = species["name"].split(" ")[0]
    if method in ("acyl_losses", "acyl_losses_nh3"):
        extra = NH3 if method == "acyl_losses_nh3" else 0.0
        hits = found(peaks_mz, precursor - FATTY_ACID_MASS - extra, tolerance)
        observed = [fa for fa, hit in zip(FATTY_ACIDS, hits) if hit]
        if not observed:
            return None
        n_chains = 3 if lipid_class == "TG" else 2
        chains = combine_chains(observed, n_chains, total)
        evidence = "fatty acid losses: " + ", ".join(chain_name(*fa) for fa in observed)
        if chains is None:  # losses seen, but no single composition: stays at species level
            return {"name": species["name"], "evidence": evidence}
        return {"name": f"{lipid_class} " + "_".join(chain_name(*chain) for chain in chains), "evidence": evidence}

    if method == "sphingoid_base":
        hits = found(peaks_mz, SPHINGOID_ION_MZ, tolerance).any(axis=1)
        options = []
        for (c, d), hit in zip(SPHINGOID_BASES, hits):
            acyl = (total[0] - c, total[1] - d)
            if hit and acyl[0] >= 2 and 0 <= acyl[1] <= acyl[0] // 2:
                options.append(((c, d), acyl))
        if len(options) != 1:  # none, or several bases: no single answer
            return None
        (c, d), acyl = options[0]
        return {"name": f"{lipid_class} {c}:{d};O2/{chain_name(*acyl)}",
                "evidence": f"sphingoid base {c}:{d};O2 ions"}
    raise ValueError(f"Unknown chain evidence '{method}' in the lipid rules")


# ---- Lipid names (Goslin) ------------------------------------------------------

_PARSER = None
_NAME_CACHE: dict[str, dict | None] = {}
GOSLIN_LEVELS = {
    "CATEGORY": "category", "CLASS": "class", "SPECIES": "species",
    "MOLECULAR_SPECIES": "molecular species", "SN_POSITION": "sn-position",
    "STRUCTURE_DEFINED": "structure defined", "FULL_STRUCTURE": "full structure",
    "COMPLETE_STRUCTURE": "complete structure",
}
# A species name written by these rules for classes Goslin does not know
# (e.g. "OL 34:1"): class, carbons, double bonds, optional oxygen suffix.
_SIMPLE_NAME = re.compile(r"(\S+) (\d+):(\d+)(;O\d)?")


def parse_lipid(name: str) -> dict | None:
    """
    Read a lipid name with Goslin (pygoslin), the reference parser of the
    lipid shorthand nomenclature. Returns the normalised ``name``, the
    ``species`` (sum composition, e.g. "PC 34:1" for "PC(16:0/18:1(9Z))"),
    the ``lipid_class`` (e.g. "PC", "PC-O" for ether lipids) and the Liebisch
    structural ``level``; None if the name is not a lipid name.

    Names Goslin does not know but that follow the species pattern of these
    rules (ornithine lipids, DGTS) are read directly, at species level.
    """
    global _PARSER
    if name in _NAME_CACHE:
        return _NAME_CACHE[name]
    if _PARSER is None:
        from pygoslin.parser.Parser import LipidParser  # 0.6 s: only when needed
        _PARSER = LipidParser()
    try:
        lipid = _PARSER.parse(name)
    except Exception:  # Goslin raises its own exception types for unknown names
        simple = _SIMPLE_NAME.fullmatch(name.strip())
        result = None if simple is None else {
            "name": name.strip(), "species": name.strip(), "lipid_class": simple.group(1), "level": "species"}
    else:
        species = lipid.get_lipid_string(type(lipid.lipid.info.level).SPECIES)
        normalised = lipid.get_lipid_string()
        level = GOSLIN_LEVELS[lipid.lipid.info.level.name]
        if species.startswith("SE 27:1/"):
            # Goslin writes cholesteryl esters as sterol esters ("SE 27:1/18:1",
            # sn-position level); one chain, like LPC: molecular species.
            species = normalised = "CE " + species.split("/", 1)[1]
            level = "molecular species"
        head, rest = species.split(" ", 1)
        result = {"name": normalised, "species": species,
                  "lipid_class": head + ("-O" if rest.startswith("O-") else ""), "level": level}
    _NAME_CACHE[name] = result
    return result


def structural_level(name: str) -> str:
    """Liebisch structural level of a lipid name (species when unknown)."""
    parsed = parse_lipid(name)
    return parsed["level"] if parsed else "species"


# ---- Annotation ---------------------------------------------------------------

def annotate_feature(
    feature_id: int, precursor: float, spectrum, species: pd.DataFrame, rules: dict, settings: LipidSettings,
) -> list[dict]:
    """
    Candidate rows for one feature. ``spectrum`` is (m/z, intensity) or None
    for a feature without MS2.
    """
    window = precursor * settings.precursor_tolerance_ppm * 1e-6
    first, last = np.searchsorted(species["mz"].to_numpy(), [precursor - window, precursor + window])
    candidates = species.iloc[first:last]
    if candidates.empty:
        return []

    peaks_mz = np.array([])
    if spectrum is not None and len(spectrum[0]):
        mz, intensity = spectrum
        peaks_mz = np.sort(mz[intensity >= settings.min_relative_intensity * intensity.max()])

    matched, mz_only = [], []
    for _, row in candidates.iterrows():
        common = {
            "feature_id": feature_id, "source": "lipid_rules", "formula": row["formula"], "adduct": row["adduct"],
            "mz_error_ppm": (precursor - row["mz"]) / row["mz"] * 1e6, "lipid_class": row["lipid_class"],
            "score_name": "fraction of rule ions found",
        }
        if len(peaks_mz):
            result = evaluate(rules[row["lipid_class"]]["adducts"][row["adduct"]], row, precursor,
                              peaks_mz, settings.fragment_tolerance_da)
            if result["matched"]:
                matched.append(common | {
                    "name": result["name"], "lipid_name": result["name"], "score": result["score"],
                    "matched_peaks": len(result["found"]), "proposed_level": "3",
                    "evidence": to_json({"found": result["found"], "missing": result["missing"]}),
                })
                continue
        mz_only.append(common | {
            "name": row["name"], "lipid_name": row["name"], "score": 0.0, "proposed_level": "5",
            "evidence": to_json({"found": [], "note": "precursor m/z only"}),
        })

    rows = matched if matched else (mz_only if settings.mz_only_candidates else [])
    # Best first: highest score, then smallest m/z error.
    rows.sort(key=lambda r: (-r["score"], abs(r["mz_error_ppm"])))
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
        row["lipid_level"] = structural_level(row["lipid_name"])
    return rows


def annotate_lipids(features: pd.DataFrame, spectra: dict, rules: dict, settings: LipidSettings) -> pd.DataFrame:
    """Candidates for every feature (``spectra``: feature id -> (m/z, intensity))."""
    species = species_table(rules)
    rows = []
    for feature in features.itertuples():
        rows.extend(annotate_feature(int(feature.feature_id), float(feature.mz),
                                     spectra.get(int(feature.feature_id)), species, rules, settings))
    table = candidate_table(rows)
    n_matched = table.loc[table["proposed_level"] == "3", "feature_id"].nunique()
    log.info("Lipid rules: %d species of %d classes; %d features with class evidence, %d with m/z matches only",
             len(species), len(rules), n_matched, table["feature_id"].nunique() - n_matched)
    return table


def run_annotate_lipids(mgf_file, features_file, rules_file, out, settings: LipidSettings) -> None:
    """File-level entry point of the ``annotate_lipids`` rule."""
    features = pd.read_parquet(features_file)
    if settings.enabled:
        table = annotate_lipids(features, read_mgf(mgf_file), load_rules(rules_file), settings)
    else:
        table = candidate_table([])
    table.to_parquet(out, index=False)
