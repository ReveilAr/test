"""
The shared format of annotation candidates.

Every annotation source writes one table with these columns, one row per
candidate identity of a feature. Any identity column may be empty (a lipid
rule gives a name and formula but no structure; a library hit may lack a
formula). Keeping one format is what lets the harmonization compare sources
and the app show them side by side.

Confidence levels follow Schymanski et al. (2014):

    1   confirmed structure: reference standard, MS2 and retention time
    2a  probable structure: library spectrum match
    2b  probable structure: diagnostic evidence (not assigned automatically in v1)
    3   tentative candidate: in-silico match, class-level evidence, lipid rules
    4   unequivocal molecular formula
    5   exact mass of interest (m/z only)
"""

import json

import numpy as np
import pandas as pd

# Column -> pandas dtype. "object" columns hold text (or lists of numbers for
# the reference spectrum of a library hit, used by the app's mirror plot).
COLUMNS = {
    "feature_id": "int64",
    "source": "object",  # e.g. "library:inhouse", "lipid_rules"
    "rank": "int64",  # 1 = the source's best candidate for this feature
    "name": "object",
    "formula": "object",
    "adduct": "object",
    "smiles": "object",
    "inchikey": "object",
    "score": "float64",  # the source's own score
    "score_name": "object",  # e.g. "cosine", "fraction of rule ions found"
    "matched_peaks": "float64",  # NaN when not applicable
    "mz_error_ppm": "float64",
    "rt_error_s": "float64",  # NaN when the reference has no retention time
    "library_kind": "object",  # experimental / in_silico / "" (not a library)
    "lipid_class": "object",  # e.g. "PC" (lipid sources, or lipid names from libraries)
    "lipid_name": "object",  # shorthand name, e.g. "PC 34:1" or "TG 16:0_18:1_18:1"
    "lipid_level": "object",  # structural level: "species", "molecular species", ...
    "proposed_level": "object",  # the source's own Schymanski level
    "evidence": "object",  # JSON text: what supports the candidate
    "reference_mz": "object",  # library spectrum (list of m/z), for mirror plots
    "reference_intensity": "object",
}

# Levels from most to least confident, and their order.
LEVELS = ["1", "2a", "2b", "3", "4", "5"]
LEVEL_RANK = {level: rank for rank, level in enumerate(LEVELS)}


def empty_table() -> pd.DataFrame:
    """A candidate table without rows (a source that found nothing)."""
    return pd.DataFrame({column: pd.Series(dtype=dtype) for column, dtype in COLUMNS.items()})


def candidate_table(rows: list[dict]) -> pd.DataFrame:
    """
    Candidate rows (dicts with any subset of COLUMNS) -> a table with every
    column, in order, with the right types. Missing text is "", missing
    numbers NaN.
    """
    if not rows:
        return empty_table()
    table = pd.DataFrame(rows)
    unknown = set(table.columns) - set(COLUMNS)
    if unknown:
        raise ValueError(f"Unknown candidate column(s): {sorted(unknown)}")
    table = table.reindex(columns=list(COLUMNS))
    for column, dtype in COLUMNS.items():
        if dtype == "object" and column not in ("reference_mz", "reference_intensity"):
            table[column] = table[column].fillna("").astype(str)
        elif dtype != "object":
            table[column] = table[column].astype(dtype)
    bad = set(table["proposed_level"]) - set(LEVELS)
    if bad:
        raise ValueError(f"Unknown confidence level(s): {sorted(bad)}")
    return table


def worst_level(*levels: str) -> str:
    """The least confident of some levels (e.g. worst_level("2a", "3") == "3")."""
    return max(levels, key=LEVEL_RANK.__getitem__)


def to_json(evidence: dict) -> str:
    """Evidence dict -> compact JSON text (numpy numbers converted)."""
    def plain(value):
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"Not JSON serializable: {value!r}")
    return json.dumps(evidence, default=plain, separators=(",", ":"))
