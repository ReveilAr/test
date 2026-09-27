"""
Spectral libraries (matchms): harmonization, then identity search.

Three steps, so that a large library is read once per project and not at
every search:

1. ``prepare_library`` reads one library file (MSP, MGF, or JSON in the GNPS
   format) and harmonizes it, in the spirit of FragHub (Dablanc et al.
   2024), with matchms's own metadata-cleaning filters (all offline):

   * metadata: field names (``Name`` / ``TITLE`` / ``compound_name``...),
     adduct notation (``M+H`` -> ``[M+H]+``), ion mode, precursor m/z,
     retention time; values stored in the wrong field (an adduct inside the
     name...) are moved to the right one;
   * structures: missing SMILES / InChI / InChIKey / formula derived from
     each other (RDKit); for spectra with a structure, annotations that
     contradict the precursor mass are repaired (wrong adduct, salt,
     molar instead of monoisotopic mass...) or, if they cannot be, the
     spectrum is removed: a spectrum whose structure does not fit its
     precursor would give confident but wrong hits;
   * spectra kept: the project's ion mode (or unknown), with a precursor
     m/z, centroided, with at least ``min_library_peaks`` fragments after
     the same cleaning as the feature spectra (``network.spectra
     .clean_peaks``: precursor region removed, intensities normalised);
   * in-silico spectra are recognised by their metadata ("in silico",
     "predicted", "LipidBlast", "CFM-ID"...), so that a mixed download
     (e.g. from MoNA) does not pass predicted spectra off as measured ones.

   What every filter removed or changed is reported (``library_cleaning.tsv``).
2. ``combine_libraries`` merges the libraries into one table and removes
   duplicates: the same spectrum (same precursor, adduct and fragments) found
   in several libraries is kept once, from the library that gives the most
   confident level (reference standards, then experimental, then in silico,
   then the order of the library list). Summary: ``library_summary.tsv``.
3. ``search_library`` compares every feature's MS2 spectrum with the library
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
from itertools import islice
from pathlib import Path

import numpy as np
import pandas as pd
from matchms import Spectrum, SpectrumProcessor
from matchms import filtering as msfilters
from matchms.filtering import default_pipelines
from matchms.importing import load_spectra
from matchms.similarity import CosineGreedy

from atlas_ms.annotation.schema import candidate_table, empty_table, to_json
from atlas_ms.config import LibrarySearchSettings
from atlas_ms.network.spectra import clean_peaks, load_feature_spectra

log = logging.getLogger(__name__)

# Columns of a prepared library (one row per spectrum, sorted by precursor_mz).
LIBRARY_COLUMNS = ["precursor_mz", "name", "adduct", "formula", "smiles", "inchikey", "rt_s",
                   "library_id", "in_silico", "mz", "intensity"]

# Words that mark a predicted spectrum in library metadata (comments,
# instrument, source...), lower case.
IN_SILICO_WORDS = ("in silico", "in-silico", "insilico", "predicted", "theoretical", "lipidblast",
                   "cfm-id", "cfmid", "virtual spectrum")
# Metadata fields not searched for those words (names and structures).
IDENTITY_FIELDS = {"compound_name", "smiles", "inchi", "inchikey", "formula", "adduct"}
BATCH = 10_000  # spectra harmonized at a time (memory stays low for large libraries)


# ---- Custom filters (matchms style: spectrum in, spectrum or None out) --------

def keep_ion_mode(spectrum_in, polarity: str, clone: bool = True):
    """Remove spectra of the other ion mode; keep those without one (common in in-house libraries)."""
    if spectrum_in is None:
        return None
    ionmode = str(spectrum_in.get("ionmode") or "").lower()
    return None if ionmode in ("positive", "negative") and ionmode != polarity else spectrum_in


def repair_structure_annotation(spectrum_in, clone: bool = True):
    """
    For spectra with a structure (SMILES or InChI) only: matchms's repair
    filters (salts, wrong adduct, molar mass given as parent mass,
    SMILES / InChI / InChIKey that disagree), then removal if the structure
    still does not fit the precursor. Spectra without a structure pass
    unchanged (the matchms filters would drop or damage them).
    """
    if spectrum_in is None or not (spectrum_in.get("smiles") or spectrum_in.get("inchi")):
        return spectrum_in
    spectrum = spectrum_in.clone() if clone else spectrum_in
    for step in REPAIR_STEPS:
        function, kwargs = step if isinstance(step, tuple) else (step, {})
        spectrum = function(spectrum, **kwargs)
        if spectrum is None:
            return None
    spectrum = msfilters.require_parent_mass_match_smiles(spectrum, mass_tolerance=0.1)
    if spectrum is not None and spectrum.get("adduct"):
        # The repairs set the parent mass from the structure when they
        # disagree: that mass must also fit the precursor m/z and adduct.
        spectrum = msfilters.require_matching_adduct_precursor_mz_parent_mass(spectrum, tolerance=0.1)
    return spectrum


def mark_in_silico(spectrum_in, clone: bool = True):
    """Record ``in_silico`` = True when the metadata calls the spectrum predicted."""
    if spectrum_in is None:
        return None
    spectrum = spectrum_in.clone() if clone else spectrum_in
    text = " ".join(str(value).lower() for key, value in spectrum.metadata.items()
                    if key not in IDENTITY_FIELDS and isinstance(value, str))
    spectrum.set("in_silico", any(word in text for word in IN_SILICO_WORDS))
    return spectrum


def _name(step) -> str:
    return (step[0] if isinstance(step, tuple) else step).__name__


# matchms's annotation repair filters, without the one that queries PubChem
# online (derive_annotation_from_compound_name).
REPAIR_STEPS = [step for step in default_pipelines.REPAIR_ANNOTATION
                if _name(step) != "derive_annotation_from_compound_name"]


def library_filters(polarity: str, precursor_window_da: float, min_peaks: int, repair: bool) -> list:
    """The harmonization pipeline (see the module docstring), in order."""
    metadata = [
        *default_pipelines.HARMONIZE_METADATA_FIELD_NAMES,
        # Not derive_formula_from_name: it takes lipid abbreviations such as
        # "POPC" for formulas and removes them from the name.
        *[step for step in default_pipelines.DERIVE_METADATA_IN_WRONG_FIELD
          if _name(step) != "derive_formula_from_name"],
        *default_pipelines.HARMONIZE_METADATA_ENTRIES,
        msfilters.require_precursor_mz,
        (keep_ion_mode, {"polarity": polarity}),
    ]
    structures = [*default_pipelines.DERIVE_MISSING_METADATA, repair_structure_annotation] if repair else []
    peaks = [
        msfilters.remove_profiled_spectra,
        (msfilters.remove_peaks_around_precursor_mz, {"mz_tolerance": precursor_window_da}),
        (msfilters.reduce_to_number_of_peaks, {"n_max": 500}),
        msfilters.normalize_intensities,
        (msfilters.require_minimum_number_of_peaks, {"n_required": min_peaks}),
    ]
    return [*metadata, *structures, *peaks, mark_in_silico]


def _text(spectrum: Spectrum, key: str) -> str:
    value = spectrum.get(key)
    return "" if value is None else str(value).strip()


def _row(spectrum: Spectrum, rt_unit: str) -> dict:
    """One harmonized spectrum as a table row."""
    rt = spectrum.get("retention_time")
    return {
        "precursor_mz": float(spectrum.get("precursor_mz")),
        "name": _text(spectrum, "compound_name"),
        "adduct": _text(spectrum, "adduct"),
        "formula": _text(spectrum, "formula"),
        "smiles": _text(spectrum, "smiles"),
        "inchikey": _text(spectrum, "inchikey"),
        "rt_s": np.nan if rt is None else float(rt) * (60.0 if rt_unit == "min" else 1.0),
        "library_id": _text(spectrum, "spectrum_id") or _text(spectrum, "spectrumid"),
        "in_silico": bool(spectrum.get("in_silico")),
        "mz": spectrum.peaks.mz.tolist(),
        "intensity": spectrum.peaks.intensities.tolist(),
    }


def prepare_library(path: str | Path, polarity: str, precursor_window_da: float, rt_unit: str = "min",
                    min_peaks: int = 2, repair: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read and harmonize one library file (see the module docstring).
    Retention times are converted to seconds (``rt_unit``: "min" or "s").
    Returns the library table, sorted by precursor m/z, and the cleaning
    report: spectra removed and changed by each filter.
    """
    processor = SpectrumProcessor(library_filters(polarity, precursor_window_da, min_peaks, repair))
    spectra = iter(load_spectra(str(path)))
    rows, report, n_read = [], None, 0
    while batch := list(islice(spectra, BATCH)):
        n_read += len(batch)
        cleaned, batch_report = processor.process_spectra(batch, progress_bar=False, create_report=True)
        rows += [_row(spectrum, rt_unit) for spectrum in cleaned]
        frame = batch_report.to_dataframe()
        report = frame if report is None else report.add(frame, fill_value=0)
        log.info("%s: %d spectra read, %d kept so far", Path(path).name, n_read, len(rows))
    report = (report if report is not None else pd.DataFrame(
        columns=["removed spectra", "changed metadata", "changed mass spectrum"])).astype(int)
    report = report.rename_axis("step").reset_index()
    report.insert(0, "spectra read", n_read)
    library = pd.DataFrame(rows, columns=LIBRARY_COLUMNS).sort_values("precursor_mz", ignore_index=True)
    log.info("%s: %d of %d spectra kept, %d marked in silico", Path(path).name, len(library), n_read,
             int(library["in_silico"].sum()))
    return library, report


def run_prepare_library(path, out, report_out, polarity: str, precursor_window_da: float, rt_unit: str,
                        min_peaks: int, repair: bool) -> None:
    """File-level entry point of the ``prepare_library`` rule."""
    library, report = prepare_library(path, polarity, precursor_window_da, rt_unit, min_peaks, repair)
    library.to_parquet(out, index=False)
    report.to_csv(report_out, sep="\t", index=False)


# ---- Combining libraries -------------------------------------------------------

def spectrum_key(row) -> tuple:
    """
    What makes two library spectra the same: precursor m/z (0.01), adduct and
    fragments (m/z to 0.01, relative intensity to 1%). Libraries copy records
    from each other (MassBank into GNPS and MoNA, for instance), and these
    copies differ at most by rounding.
    """
    peaks = tuple(zip(np.round(row.mz, 2), np.round(row.intensity, 2)))
    return round(row.precursor_mz, 2), row.adduct, peaks


def combine_libraries(tables: list[pd.DataFrame], entries: list[dict], remove_duplicates: bool = True):
    """
    One table from all libraries, each spectrum with its ``library`` and its
    ``kind`` (in silico if its library is, or if its metadata says so), and
    duplicates removed (see the module docstring). Returns the table, sorted
    by precursor m/z, and the summary per library.
    """
    def priority(entry):  # lower first: the library whose copy is kept
        return (0 if entry["reference_standards"] else 1) + (2 if entry["kind"] == "in_silico" else 0)

    order = sorted(range(len(entries)), key=lambda i: (priority(entries[i]), i))
    parts, summary, seen = [], [], set()
    for i in order:
        table, entry = tables[i], entries[i]
        kind = np.where((entry["kind"] == "in_silico") | table["in_silico"].to_numpy(dtype=bool),
                        "in_silico", "experimental")
        table = table.assign(library=entry["name"], kind=kind)
        keep = np.ones(len(table), dtype=bool)
        if remove_duplicates:
            for position, row in enumerate(table.itertuples()):
                key = spectrum_key(row)
                keep[position] = key not in seen
                seen.add(key)
        parts.append(table[keep])
        summary.append({"library": entry["name"], "spectra after cleaning": len(table),
                        "marked in silico": int((kind == "in_silico").sum()) if entry["kind"] != "in_silico" else 0,
                        "duplicates removed": int((~keep).sum()), "spectra searched": int(keep.sum())})
    combined = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
        {column: [] for column in [*LIBRARY_COLUMNS, "library", "kind"]})
    summary = pd.DataFrame(summary, columns=["library", "spectra after cleaning", "marked in silico",
                                             "duplicates removed", "spectra searched"])
    # Back to the order of the library list, for reading.
    summary = summary.set_index("library").loc[[e["name"] for e in entries]].reset_index()
    return combined.sort_values("precursor_mz", ignore_index=True), summary


def run_combine_libraries(library_files: list, report_files: list, entries: list[dict], out,
                          summary_out, cleaning_out, remove_duplicates: bool) -> None:
    """File-level entry point of the ``combine_libraries`` rule."""
    tables = [pd.read_parquet(path) for path in library_files]
    combined, summary = combine_libraries(tables, entries, remove_duplicates)
    reports = [pd.read_csv(path, sep="\t").assign(library=entry["name"]) for path, entry in zip(report_files, entries)]
    cleaning = pd.concat(reports, ignore_index=True) if reports else pd.DataFrame(columns=["library", "spectra read"])
    read = cleaning.groupby("library")["spectra read"].first()
    summary.insert(1, "spectra read", summary["library"].map(read))
    combined.to_parquet(out, index=False)
    summary.to_csv(summary_out, sep="\t", index=False)
    cleaning.to_csv(cleaning_out, sep="\t", index=False)
    log.info("Libraries combined: %d spectra searched (%d duplicates removed)",
             len(combined), int(summary["duplicates removed"].sum()))


# ---- Search ------------------------------------------------------------------

def search_library(
    queries: list[Spectrum],
    feature_rt: dict[int, float],
    library: pd.DataFrame,
    entries: dict[str, dict],
    settings: LibrarySearchSettings,
) -> pd.DataFrame:
    """
    Identity search of the feature spectra (``queries``, cleaned, with a
    ``feature_id``) in the combined library (``combine_libraries``: each
    spectrum with its ``library`` and ``kind``). ``entries`` gives each
    library's settings by name (``LibrarySearchSettings.library``). Returns
    candidates in the shared format, at most ``top_n`` per feature and
    library (so that a large public library cannot crowd out the hit of a
    small library of standards).
    """
    cosine = CosineGreedy(tolerance=settings.fragment_tolerance_da)
    # Plain arrays: much faster than reading the table row by row.
    precursors = library["precursor_mz"].to_numpy()
    peak_mz, peak_intensity = library["mz"].to_numpy(), library["intensity"].to_numpy()
    names = library["library"].to_numpy()
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
        for name in dict.fromkeys(names[index] for _, _, index in hits):  # each library, best hit first
            best = [hit for hit in hits if names[hit[2]] == name][:settings.top_n]
            for rank, (score, matches, index) in enumerate(best, start=1):
                rows.append(_candidate(feature_id, precursor, feature_rt.get(feature_id, np.nan),
                                       rank, score, matches, library.iloc[index], entries[name]))
    return candidate_table(rows)


def _candidate(feature_id, precursor, rt, rank, score, matches, reference, entry) -> dict:
    """One library hit as a candidate row, with its proposed confidence level."""
    rt_error = rt - reference["rt_s"]  # NaN when the library has no retention time
    rt_ok = bool(np.abs(rt_error) <= entry["rt_tolerance_s"])  # False for NaN
    if reference["kind"] == "in_silico":  # the library's kind, or the spectrum's own metadata
        level = "3"
    elif entry["reference_standards"] and rt_ok:
        level = "1"
    else:
        level = "2a"
    evidence = {"library_id": reference["library_id"], "reference_precursor_mz": reference["precursor_mz"]}
    if entry["reference_standards"]:
        evidence["retention_time"] = "matches" if rt_ok else "does not match or unknown"
    if reference["kind"] == "in_silico" and entry["kind"] != "in_silico":
        evidence["in_silico"] = "marked as predicted in the library's metadata"
    return {
        "feature_id": feature_id, "source": f"library:{entry['name']}", "rank": rank,
        "name": reference["name"], "formula": reference["formula"], "adduct": reference["adduct"],
        "smiles": reference["smiles"], "inchikey": reference["inchikey"],
        "score": score, "score_name": "cosine", "matched_peaks": matches,
        "mz_error_ppm": (precursor - reference["precursor_mz"]) / reference["precursor_mz"] * 1e6,
        "rt_error_s": rt_error, "library_kind": reference["kind"], "proposed_level": level,
        "evidence": to_json(evidence),
        "reference_mz": list(reference["mz"]), "reference_intensity": list(reference["intensity"]),
    }


def run_search_libraries(
    mgf_file, features_file, library_file, entries: list[dict], out,
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
    library = pd.read_parquet(library_file) if library_file else pd.DataFrame()
    if library.empty:
        result = empty_table()
    else:
        result = search_library(queries, feature_rt, library, {e["name"]: e for e in entries}, settings)
    log.info("Library search: %d hits for %d of %d spectra, from %s", len(result), result["feature_id"].nunique(),
             len(queries), result["source"].value_counts().to_dict())
    result.to_parquet(out, index=False)
