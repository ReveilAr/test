"""
MS2Query (de Jonge et al. 2023), in its own conda environment.

MS2Query searches a large library of reference spectra (GNPS-derived, in
its own files) for the best *analog* of each spectrum: the same compound or
a close structural relative. Its score (the "metascore", 0-1) combines
MS2DeepScore, Spec2Vec and the precursor mass difference; 0.7 and above is
where the authors report reliable analogs.

It needs matchms <= 0.26.4, ms2deepscore 2.0.0 and torch < 2.6, which
conflict with the core environment: this module runs in the ``ms2query``
environment (``workflow/envs/ms2query.yaml``) and imports ms2query only
inside the functions that use it.

Candidates (``source`` "ms2query") propose level 3. An analog with the same
precursor m/z is probably the compound itself, but MS2Query reports no
cosine score or matched fragments to confirm it (the level 2a rule), so
exact matches also stay at level 3; the evidence says they are exact.
"""

import json
import logging
from pathlib import Path

import pandas as pd

from atlas_ms.annotation.schema import candidate_table, to_json

log = logging.getLogger(__name__)

# The library and models (positive mode, a few GB), downloaded once per machine.
DEFAULT_MODELS_DIR = Path.home() / ".cache" / "atlas-ms" / "ms2query" / "positive"
DOWNLOAD_DONE = DEFAULT_MODELS_DIR / "download.done"  # written when every file is there
CLASS_COLUMNS = ["cf_kingdom", "cf_superclass", "cf_class", "cf_subclass", "cf_direct_parent"]


def models_input(models_dir: str) -> str:
    """What the MS2Query rule waits for: your own models folder, or the downloaded one."""
    return models_dir or str(DOWNLOAD_DONE)


def parse_results(csv_file, precursor_tolerance_ppm: float) -> pd.DataFrame:
    """
    MS2Query's CSV -> candidates. The feature id comes from the MGF
    ``SCANS`` field (asked for as an extra column). Rows of one feature are
    already sorted by score, best first.
    """
    table = pd.read_csv(csv_file)
    rows = []
    for feature_id, group in table.groupby("scans", sort=False):
        for rank, hit in enumerate(group.itertuples(), start=1):
            precursor = float(hit.precursor_mz_query_spectrum)
            exact = abs(hit.precursor_mz_difference) <= precursor * precursor_tolerance_ppm * 1e-6
            classes = [str(getattr(hit, column)) for column in CLASS_COLUMNS
                       if column in table and isinstance(getattr(hit, column), str)]
            evidence = {"exact_precursor_match": bool(exact),
                        "precursor_mz_difference": round(float(hit.precursor_mz_difference), 4)}
            if classes:
                evidence["classyfire"] = " > ".join(classes)
            if "npc_class_results" in table and isinstance(hit.npc_class_results, str):
                evidence["npc_class"] = hit.npc_class_results
            rows.append({
                "feature_id": int(feature_id), "source": "ms2query", "rank": rank,
                "name": hit.analog_compound_name if isinstance(hit.analog_compound_name, str) else hit.inchikey,
                "smiles": hit.smiles if isinstance(hit.smiles, str) else "",
                "inchikey": hit.inchikey if isinstance(hit.inchikey, str) else "",
                "score": float(hit.ms2query_model_prediction), "score_name": "MS2Query score",
                "mz_error_ppm": (precursor - hit.precursor_mz_analog) / hit.precursor_mz_analog * 1e6,
                "proposed_level": "3", "evidence": to_json(evidence),
            })
    return candidate_table(rows)


def run_ms2query(mgf_file, models_dir, csv_out, out, settings: dict, threads: int = 1) -> None:
    """
    File-level entry point of the ``run_ms2query`` rule. MS2Query writes its
    CSV into ``csv_out``'s folder, named after the MGF file (and adds a
    number to the name if that file exists: it is deleted first).
    """
    import torch
    from ms2query.ms2library import create_library_object_from_one_dir
    from ms2query.run_ms2query import run_ms2query_single_file
    from ms2query.utils import SettingsRunMS2Query

    torch.set_num_threads(threads)  # MS2DeepScore inside MS2Query: the cores Snakemake gave this step

    models_dir = Path(models_dir)
    if models_dir.is_file():  # the download marker: the files are next to it
        models_dir = models_dir.parent
    library = create_library_object_from_one_dir(str(models_dir))
    csv_out = Path(csv_out)
    written = csv_out.parent / (Path(mgf_file).stem + ".csv")
    written.unlink(missing_ok=True)
    run_ms2query_single_file(
        spectrum_file_name=Path(mgf_file).name, folder_with_spectra=str(Path(mgf_file).parent),
        results_folder=str(csv_out.parent), ms2library=library,
        settings=SettingsRunMS2Query(
            nr_of_top_analogs_to_save=settings["top_n"], minimal_ms2query_metascore=settings["min_score"],
            additional_metadata_columns=("scans",), filter_on_ion_mode="positive"),
    )
    written = csv_out.parent / (Path(mgf_file).stem + ".csv")
    if written != csv_out:
        written.rename(csv_out)
    if csv_out.exists():
        candidates = parse_results(csv_out, settings["precursor_tolerance_ppm"])
    else:  # no spectrum had an analog above the minimum score
        candidates = candidate_table([])
    log.info("MS2Query: %d candidates for %d features", len(candidates), candidates["feature_id"].nunique())
    candidates.to_parquet(out, index=False)


def download_models(done_file) -> None:
    """Download MS2Query's positive-mode library and models (once per machine)."""
    from ms2query.run_ms2query import download_zenodo_files

    directory = Path(done_file).parent
    directory.mkdir(parents=True, exist_ok=True)
    download_zenodo_files("positive", str(directory))
    Path(done_file).write_text(json.dumps(sorted(p.name for p in directory.iterdir())))
