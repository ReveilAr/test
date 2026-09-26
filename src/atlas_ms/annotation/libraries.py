"""
Spectral library search (matchms).

Two steps, so that a large library is read once per project and not at
every search:

1. ``prepare_library`` reads one library file (MSP, MGF, or JSON in the GNPS
   format) with matchms, keeps the spectra of the project's ion mode, cleans
   them exactly like the feature spectra (``network.spectra.clean_peaks``:
   precursor region removed, intensities normalised) and stores them as a
   Parquet table sorted by precursor m/z.
2. ``search_library`` compares every feature's MS2 spectrum with the library
   spectra whose precursor m/z is within ``precursor_tolerance_ppm`` (an
   *identity* search: same precursor, similar fragments). The score is the
   cosine (matchms ``CosineGreedy``), which also counts matched fragments.

Confidence levels proposed for a hit (see ``schema.py``):

* 1: library of reference standards measured on the same method, and the
  retention time matches (``rt_tolerance_s``);
* 2a: any other experimental spectrum (including a reference standard whose
  retention time is off or unknown);
* 3: in-silico (predicted) spectra, e.g. LipidBlast: a good match with a
  predicted spectrum is weaker evidence than with a measured one.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from matchms import Spectrum
from matchms import filtering as msfilters
from matchms.importing import load_spectra
from matchms.similarity import CosineGreedy

from atlas_ms.annotation.schema import candidate_table, empty_table, to_json
from atlas_ms.config import LibrarySearchSettings
from atlas_ms.network.spectra import clean_peaks, load_feature_spectra

log = logging.getLogger(__name__)

# Columns of a prepared library (one row per spectrum, sorted by precursor_mz).
LIBRARY_COLUMNS = ["precursor_mz", "name", "adduct", "formula", "smiles", "inchikey", "rt_s",
                   "library_id", "mz", "intensity"]


def _text(spectrum: Spectrum, key: str) -> str:
    value = spectrum.get(key)
    return "" if value is None else str(value).strip()


def prepare_library(path: str | Path, polarity: str, precursor_window_da: float, rt_unit: str = "min") -> pd.DataFrame:
    """
    Read and clean one library file. Spectra of the other ion mode, without
    precursor m/z or without fragments left after cleaning are skipped.
    Spectra without an ion mode are kept (many in-house libraries omit it).
    Retention times are converted to seconds (``rt_unit``: "min" or "s").
    """
    rows, skipped = [], {"ion mode": 0, "no precursor m/z": 0, "no fragments": 0}
    for spectrum in load_spectra(str(path)):
        spectrum = msfilters.default_filters(spectrum)  # harmonised metadata (name, adduct, precursor m/z...)
        ionmode = _text(spectrum, "ionmode").lower()
        if ionmode in ("positive", "negative") and ionmode != polarity:
            skipped["ion mode"] += 1
            continue
        spectrum = msfilters.require_precursor_mz(spectrum)
        if spectrum is None:
            skipped["no precursor m/z"] += 1
            continue
        spectrum = msfilters.add_retention_time(spectrum)
        cleaned = clean_peaks(spectrum, precursor_window_da)
        if cleaned is None or len(cleaned.peaks.mz) == 0:
            skipped["no fragments"] += 1
            continue
        rt = spectrum.get("retention_time")
        rows.append({
            "precursor_mz": float(spectrum.get("precursor_mz")),
            "name": _text(spectrum, "compound_name"),
            "adduct": _text(spectrum, "adduct"),
            "formula": _text(spectrum, "formula"),
            "smiles": _text(spectrum, "smiles"),
            "inchikey": _text(spectrum, "inchikey"),
            "rt_s": np.nan if rt is None else float(rt) * (60.0 if rt_unit == "min" else 1.0),
            "library_id": _text(spectrum, "spectrum_id") or _text(spectrum, "spectrumid"),
            "mz": cleaned.peaks.mz.tolist(),
            "intensity": cleaned.peaks.intensities.tolist(),
        })
    log.info("%s: %d spectra kept; skipped: %s", Path(path).name, len(rows), skipped)
    library = pd.DataFrame(rows, columns=LIBRARY_COLUMNS)
    return library.sort_values("precursor_mz", ignore_index=True)


def run_prepare_library(path, out, polarity: str, precursor_window_da: float, rt_unit: str) -> None:
    """File-level entry point of the ``prepare_library`` rule."""
    prepare_library(path, polarity, precursor_window_da, rt_unit).to_parquet(out, index=False)


def search_library(
    queries: list[Spectrum],
    feature_rt: dict[int, float],
    library: pd.DataFrame,
    entry: dict,
    settings: LibrarySearchSettings,
) -> pd.DataFrame:
    """
    Identity search of the feature spectra (``queries``, cleaned, with a
    ``feature_id``) in one prepared library. ``entry`` is the library's
    settings (``LibrarySearchSettings.library``). Returns candidates in the
    shared format, at most ``top_n`` per feature.
    """
    cosine = CosineGreedy(tolerance=settings.fragment_tolerance_da)
    # Plain arrays: much faster than reading the table row by row.
    precursors = library["precursor_mz"].to_numpy()
    peak_mz, peak_intensity = library["mz"].to_numpy(), library["intensity"].to_numpy()
    rows = []
    for query in queries:
        feature_id = int(query.get("feature_id"))
        precursor = float(query.get("precursor_mz"))
        window = precursor * settings.precursor_tolerance_ppm * 1e-6
        first, last = np.searchsorted(precursors, [precursor - window, precursor + window])
        hits = []
        for index in range(first, last):  # the library is sorted: only spectra in the window
            spectrum = Spectrum(np.asarray(peak_mz[index], dtype=float), np.asarray(peak_intensity[index], dtype=float),
                                metadata={"precursor_mz": precursors[index]}, metadata_harmonization=False)
            result = cosine.pair(query, spectrum)
            score, matches = float(result["score"]), int(result["matches"])
            if score >= settings.min_score and matches >= settings.min_matched_peaks:
                hits.append((score, matches, index))
        hits.sort(key=lambda hit: -hit[0])
        for rank, (score, matches, index) in enumerate(hits[:settings.top_n], start=1):
            rows.append(_candidate(feature_id, precursor, feature_rt.get(feature_id, np.nan),
                                   rank, score, matches, library.iloc[index], entry))
    return candidate_table(rows)


def _candidate(feature_id, precursor, rt, rank, score, matches, reference, entry) -> dict:
    """One library hit as a candidate row, with its proposed confidence level."""
    rt_error = rt - reference["rt_s"]  # NaN when the library has no retention time
    rt_ok = bool(np.abs(rt_error) <= entry["rt_tolerance_s"])  # False for NaN
    if entry["kind"] == "in_silico":
        level = "3"
    elif entry["reference_standards"] and rt_ok:
        level = "1"
    else:
        level = "2a"
    evidence = {"library_id": reference["library_id"], "reference_precursor_mz": reference["precursor_mz"]}
    if entry["reference_standards"]:
        evidence["retention_time"] = "matches" if rt_ok else "does not match or unknown"
    return {
        "feature_id": feature_id, "source": f"library:{entry['name']}", "rank": rank,
        "name": reference["name"], "formula": reference["formula"], "adduct": reference["adduct"],
        "smiles": reference["smiles"], "inchikey": reference["inchikey"],
        "score": score, "score_name": "cosine", "matched_peaks": matches,
        "mz_error_ppm": (precursor - reference["precursor_mz"]) / reference["precursor_mz"] * 1e6,
        "rt_error_s": rt_error, "library_kind": entry["kind"], "proposed_level": level,
        "evidence": to_json(evidence),
        "reference_mz": list(reference["mz"]), "reference_intensity": list(reference["intensity"]),
    }


def run_search_libraries(
    mgf_file, features_file, library_files: list, entries: list[dict], out,
    precursor_window_da: float, settings: LibrarySearchSettings,
) -> None:
    """
    File-level entry point of the ``search_libraries`` rule. Every feature
    spectrum with at least one fragment is searched (no minimum peak count:
    a sparse lipid spectrum can still match its library spectrum).
    """
    queries = []
    for spectrum in load_feature_spectra(mgf_file):
        spectrum = msfilters.default_filters(spectrum)
        cleaned = clean_peaks(spectrum, precursor_window_da)
        if cleaned is not None and len(cleaned.peaks.mz):
            queries.append(cleaned)
    feature_rt = pd.read_parquet(features_file).set_index("feature_id")["rt"].to_dict()

    tables = []
    for path, entry in zip(library_files, entries):
        found = search_library(queries, feature_rt, pd.read_parquet(path), entry, settings)
        log.info("Library %s: %d hits for %d of %d spectra", entry["name"], len(found),
                 found["feature_id"].nunique(), len(queries))
        tables.append(found)
    tables = [table for table in tables if not table.empty]
    result = pd.concat(tables, ignore_index=True) if tables else empty_table()
    result.to_parquet(out, index=False)
