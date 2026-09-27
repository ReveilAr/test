"""
SIRIUS 6 through its REST API (PySirius), in its own conda environment.

This module runs in the ``sirius`` environment (``workflow/envs/sirius.yaml``),
which has SIRIUS and its Python client but not pyOpenMS: it only imports the
candidate format (``schema.py``) and gets PySirius passed in, so that the
tests can replace SIRIUS with a fake one.

One SIRIUS job covers, for every feature of ``work/sirius/input.json``
(``sirius_input.py``):

* molecular formula identification (SIRIUS: isotope pattern +
  fragmentation tree), re-ranked by ZODIAC with the whole data set;
* El Gordo: recognises lipid spectra and names the lipid species;
* CSI:FingerID: predicts the molecular fingerprint and ranks the structures
  of the chosen databases, with the COSMIC confidence of the top hit, and
  (expansive search) falls back to PubChem when the confidence is low;
* MSNovelist (optional): generates structures de novo from the predicted
  fingerprint, for compounds missing from every database;
* CANOPUS: compound classes (ClassyFire, NPClassifier).

Candidates written, and the levels they propose:

* ``sirius:formula``: formula candidates; the top one is level 4 when its
  ZODIAC score (SIRIUS score without ZODIAC) reaches ``min_formula_score``,
  the others 5;
* ``sirius:elgordo``: lipid species, level 3;
* ``sirius:csi``: structure candidates, level 3 whatever the COSMIC
  confidence (a design decision: CSI:FingerID is never more than tentative);
* ``sirius:msnovelist``: de novo structures, level 3;
* ``sirius:canopus``: the most specific compound class, level 3.

The SIRIUS project (``work/sirius/project.sirius``) is kept: it can be opened
in the SIRIUS GUI to look at fragmentation trees and every candidate.

Login: CSI:FingerID, CANOPUS and MSNovelist are web services of Bright Giant
and need a (free, academic) account. When SIRIUS is not logged in, the
account saved in the app (``atlas_ms.credentials``) is used.
"""

import hashlib
import json
import logging
import signal
import time
from pathlib import Path

from atlas_ms.annotation.schema import candidate_table, to_json
from atlas_ms.credentials import load_sirius_account

log = logging.getLogger(__name__)

PROFILES = {"orbitrap": "ORBITRAP", "qtof": "QTOF"}
FINISHED = {"DONE", "FAILED", "CANCELED"}


# ---- SIRIUS input -------------------------------------------------------------

def feature_import(models, entry: dict):
    """One feature of input.json as a SIRIUS FeatureImport (our feature id as external id)."""
    def spectrum(peaks, ms_level, name):
        return models.BasicSpectrum(
            name=name, ms_level=ms_level, precursor_mz=entry["mz"] if ms_level == 2 else None, cosine_query=False,
            peaks=[models.SimplePeak(mz=mz, intensity=intensity) for mz, intensity in peaks],
        )
    return models.FeatureImport(
        name=f"feature {entry['feature_id']}",
        external_feature_id=str(entry["feature_id"]),
        ion_mass=entry["mz"],
        charge=entry["charge"],
        detected_adducts=entry["adducts"] or None,
        rt_apex_seconds=entry["rt_s"],
        merged_ms1=spectrum(entry["ms1"], 1, "isotope pattern") if entry["ms1"] else None,
        ms2_spectra=[spectrum(entry["ms2"], 2, "MS2")],
    )


def job_submission(models, settings: dict, fallback_adducts: list[str]):
    """The SIRIUS job: formula + ZODIAC, El Gordo, CSI:FingerID, CANOPUS (see the module docstring)."""
    return models.JobSubmission(
        fallback_adducts=fallback_adducts,  # tried for features whose adduct is unknown
        formula_id_params=models.Sirius(
            enabled=True, profile=PROFILES[settings["profile"]],
            number_of_candidates=settings["formula_candidates"],
            mass_accuracy_ms2ppm=settings["ms2_ppm"],
            enforce_el_gordo_formula=True,  # lipids: use El Gordo's formula
        ),
        zodiac_params=models.Zodiac(enabled=settings["zodiac"]),
        fingerprint_prediction_params=models.FingerprintPrediction(enabled=True),
        canopus_params=models.Canopus(enabled=settings["canopus"]),
        structure_db_search_params=models.StructureDbSearch(
            enabled=True, structure_search_dbs=settings["structure_databases"],
            tag_structures_with_lipid_class=True,
            expansive_search_confidence_mode=settings["expansive_search"],
        ),
        ms_novelist_params=models.MsNovelist(
            enabled=settings["msnovelist"], number_of_candidate_to_predict=settings["msnovelist_candidates"]),
        spectra_search_params=models.SpectralLibrarySearch(enabled=False),  # our own library search does it
    )


def project_id(project_file: Path) -> str:
    """A SIRIUS project id ([a-zA-Z0-9_-]) unique to this ATLAS-MS project."""
    return "atlas_ms_" + hashlib.md5(str(Path(project_file).resolve()).encode()).hexdigest()[:10]


def wait_for(api, pid: str, job, poll_s: float = 5.0) -> None:
    """Poll the job until it finishes; log its progress every minute; raise if it failed."""
    last_log = 0.0
    while True:
        progress = api.jobs().get_job(pid, job.id).progress
        if progress.state in FINISHED:
            break
        if time.monotonic() - last_log > 60:
            log.info("SIRIUS: %s %s/%s %s", progress.state, progress.current_progress, progress.max_progress,
                     progress.message or "")
            last_log = time.monotonic()
        time.sleep(poll_s)
    if progress.state != "DONE":
        raise RuntimeError(f"SIRIUS job {progress.state}: {progress.error_message}")


# ---- SIRIUS results ------------------------------------------------------------

def formula_rows(feature_id: int, formulas: list, settings: dict) -> list[dict]:
    """Formula candidates (sirius:formula) and El Gordo lipid species (sirius:elgordo)."""
    rows = []
    for rank, formula in enumerate(formulas[:settings["candidates_kept"]], start=1):
        score, score_name = formula.zodiac_score, "ZODIAC score"
        if score is None:
            score, score_name = formula.sirius_score_normalized, "SIRIUS score (normalized)"
        confident = rank == 1 and score is not None and score >= settings["min_formula_score"]
        common = {"feature_id": feature_id, "rank": rank, "formula": formula.molecular_formula,
                  "adduct": formula.adduct, "score": score, "score_name": score_name}
        rows.append(common | {
            "source": "sirius:formula", "name": formula.molecular_formula, "proposed_level": "4" if confident else "5",
            "evidence": to_json({"explained_peaks": f"{formula.num_of_explained_peaks}/{formula.num_of_explainable_peaks}",
                                 "sirius_score": formula.sirius_score}),
        })
        lipid = formula.lipid_annotation
        if lipid is not None and lipid.lipid_species:
            rows.append(common | {
                "source": "sirius:elgordo", "name": lipid.lipid_species, "proposed_level": "3",
                "evidence": to_json({"lipid_class": lipid.lipid_class_name, "lipid_maps_id": lipid.lipid_maps_id,
                                     "hypothetical_structure": lipid.hypothetical_structure}),
            })
    return rows


def structure_rows(feature_id: int, structures: list, annotations, settings: dict,
                   source: str = "sirius:csi") -> list[dict]:
    """
    Structure candidates: CSI:FingerID's (``sirius:csi``), with the COSMIC
    confidence of the top hit, or MSNovelist's (``sirius:msnovelist``).
    """
    rows = []
    for rank, structure in enumerate(structures[:settings["candidates_kept"]], start=1):
        evidence = {"tanimoto": structure.tanimoto_similarity,
                    "databases": sorted({link.name for link in (structure.db_links or []) if link.name})}
        if rank == 1 and annotations is not None and source == "sirius:csi":
            evidence["cosmic_confidence_exact"] = annotations.confidence_exact_match
            evidence["cosmic_confidence_approximate"] = annotations.confidence_approx_match
        rows.append({
            "feature_id": feature_id, "source": source, "rank": rank,
            "name": structure.structure_name or structure.inchi_key, "smiles": structure.smiles,
            "inchikey": structure.inchi_key, "formula": structure.molecular_formula, "adduct": structure.adduct,
            "score": structure.csi_score, "score_name": "CSI:FingerID score", "proposed_level": "3",
            "evidence": to_json(evidence),
        })
    return rows


def class_rows(feature_id: int, classes) -> list[dict]:
    """The most specific CANOPUS class (sirius:canopus), with the whole lineage as evidence."""
    lineage = [c for c in (classes.classy_fire_lineage or []) if c is not None]
    best = lineage[-1] if lineage else classes.npc_class
    if best is None:
        return []
    evidence = {"classyfire": " > ".join(c.name for c in lineage)}
    for key, value in (("npc_pathway", classes.npc_pathway), ("npc_superclass", classes.npc_superclass),
                       ("npc_class", classes.npc_class)):
        if value is not None:
            evidence[key] = value.name
    return [{"feature_id": feature_id, "source": "sirius:canopus", "rank": 1, "name": best.name,
             "score": best.probability, "score_name": "CANOPUS probability", "proposed_level": "3",
             "evidence": to_json(evidence)}]


def collect_candidates(api, pid: str, settings: dict):
    """Every feature's SIRIUS results, as candidates (in the shared format)."""
    features_api = api.features()
    rows = []
    for aligned in features_api.get_aligned_features(pid, opt_fields=["topAnnotations"]):
        feature_id, afid = int(aligned.external_feature_id), aligned.aligned_feature_id
        formulas = features_api.get_formula_candidates(pid, afid, opt_fields=["statistics", "lipidAnnotation"])
        rows += formula_rows(feature_id, formulas, settings)
        structures = features_api.get_structure_candidates(pid, afid, opt_fields=["dbLinks"])
        rows += structure_rows(feature_id, structures, aligned.top_annotations, settings)
        if settings["msnovelist"]:
            de_novo = features_api.get_de_novo_structure_candidates(pid, afid, opt_fields=["dbLinks"])
            rows += structure_rows(feature_id, de_novo, None, settings, source="sirius:msnovelist")
        if settings["canopus"] and formulas:
            try:
                rows += class_rows(feature_id, features_api.get_best_matching_compound_classes(
                    pid, afid, formulas[0].formula_id))
            except Exception as error:  # no CANOPUS result for this feature (e.g. no fingerprint)
                log.debug("No compound classes for feature %s: %s", feature_id, error)
    return candidate_table(rows)


# ---- The whole run ---------------------------------------------------------------

def ensure_login(api, models, account: dict | None) -> None:
    """
    Log SIRIUS in with the saved account if it is not logged in already.
    SIRIUS keeps its session afterwards, so this usually happens once.
    """
    if api.account().is_logged_in():
        return
    if account is None:
        raise RuntimeError(
            "SIRIUS is not logged in. CSI:FingerID, CANOPUS and MSNovelist need a (free, academic) SIRIUS "
            "account: enter it in the app (Setup tab, SIRIUS card), or log in once in the SIRIUS GUI.")
    if not account.get("accept_terms"):
        raise RuntimeError("Accept the SIRIUS terms of service in the app (Setup tab, SIRIUS card) to log in.")
    try:
        api.account().login(accept_terms=True, account_credentials=models.AccountCredentials(
            username=account["username"], password=account["password"]))
    except Exception as error:  # PySirius raises an ApiException with the server's message
        raise RuntimeError(f"SIRIUS login failed for {account['username']}: {error}") from None
    if not api.account().is_logged_in():
        raise RuntimeError(f"SIRIUS login failed for {account['username']}: check the account in the Setup tab.")
    log.info("SIRIUS: logged in as %s", account["username"])


def annotate_with_sirius(api, models, entries: list[dict], settings: dict, fallback_adducts: list[str],
                         project_file: Path, poll_s: float = 5.0, account: dict | None = None):
    """
    Import the features into a new SIRIUS project, run the job, and return
    the candidates. ``api`` is a connected PySirius API, ``models`` the
    PySirius module (or fakes, in tests), ``account`` the saved SIRIUS
    account (used only if SIRIUS is not logged in).
    """
    ensure_login(api, models, account)
    pid = project_id(project_file)
    projects = api.projects()
    if any(p.project_id == pid for p in projects.get_projects()):
        projects.close_project(pid)  # left open by an earlier run
    project_file = Path(project_file)
    if project_file.exists():
        project_file.unlink()  # a re-run starts from scratch
    projects.create_project(pid, path_to_project=str(project_file.resolve()))
    try:
        profile = PROFILES[settings["profile"]]
        api.features().add_aligned_features(pid, [feature_import(models, e) for e in entries], profile=profile)
        log.info("SIRIUS: %d features imported, job started", len(entries))
        job = api.jobs().start_job(pid, job_submission(models, settings, fallback_adducts))
        wait_for(api, pid, job, poll_s)
        candidates = collect_candidates(api, pid, settings)
    finally:
        projects.close_project(pid)
    log.info("SIRIUS: %d candidates for %d features", len(candidates), candidates["feature_id"].nunique())
    return candidates


def _stop(signum, frame):
    """SIGTERM handler: leave through the normal exit path (runs the `finally` blocks)."""
    raise SystemExit(f"Stopped (signal {signum})")


def run_sirius_job(input_file, out, project_file, settings: dict, fallback_adducts: list[str]) -> None:
    """
    File-level entry point of the ``run_sirius`` rule: attach to a running
    SIRIUS 6 (e.g. your GUI) or start a headless one, run, write the
    candidates, and stop SIRIUS again if it was started here.
    """
    import PySirius  # only in the sirius environment
    from PySirius import SiriusSDK

    entries = json.loads(Path(input_file).read_text())
    sdk = SiriusSDK()
    api = sdk.attach_or_start_sirius(headless=True)
    if api is None:
        raise RuntimeError("Could not start or attach to SIRIUS 6 (see the messages above).")
    # Stop (the app's Stop button, Ctrl+C) sends SIGTERM: turn it into an
    # exception, so that the `finally` below still shuts SIRIUS down instead
    # of leaving it running in the background.
    signal.signal(signal.SIGTERM, _stop)
    try:
        candidates = annotate_with_sirius(api, PySirius, entries, settings, fallback_adducts, Path(project_file),
                                          account=load_sirius_account())
    finally:
        if SiriusSDK.process is not None:  # started by this rule, not the user's own SIRIUS
            sdk.shutdown_sirius()
    candidates.to_parquet(out, index=False)
