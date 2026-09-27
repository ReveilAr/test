"""
SIRIUS and MS2Query. Neither runs here (no SIRIUS account, no network access
for MS2Query's models): SIRIUS is replaced by a fake client with the same
methods as PySirius, MS2Query by a results file in its format. The real
tools run on your machine.
"""

import json
import shutil
from types import SimpleNamespace

import pandas as pd
import pytest
import yaml

from atlas_ms.annotation.ms2query import parse_results
from atlas_ms.annotation.sirius import annotate_with_sirius, wait_for
from atlas_ms.annotation.sirius_input import sirius_input
from atlas_ms.mgf import read_mgf
from atlas_ms.preprocessing.msdata import RunReader
from conftest import snakemake
from test_workflow import feature_of, planned_jobs

SETTINGS = {"profile": "orbitrap", "ms2_ppm": 5.0, "formula_candidates": 10, "zodiac": True,
            "min_formula_score": 0.9, "structure_databases": ["BIO"], "expansive_search": "APPROXIMATE",
            "canopus": True, "candidates_kept": 3, "msnovelist": True, "msnovelist_candidates": 64}


# ---- SIRIUS input ------------------------------------------------------------

def test_sirius_input_has_isotope_patterns_and_adducts(processed_project):
    root = processed_project.root
    features = pd.read_parquet(root / "results" / "features.parquet")
    spectra = read_mgf(root / "results" / "gnps" / "ms2_spectra.mgf")
    entries = sirius_input(features, pd.read_parquet(root / "results" / "quant.parquet"), spectra,
                           RunReader(root, list(processed_project.load_samples().index)), "positive", 10.0)
    assert {e["feature_id"] for e in entries} == set(spectra)  # every feature with MS2
    pc = next(e for e in entries if e["feature_id"] == feature_of(features, "PC 34:1 [M+H]+")["feature_id"])
    assert pc["charge"] == 1 and pc["adducts"] == ["[M+H]+"]
    # The synthetic isotope pattern: M, M+1 (45%), M+2 (12%).
    assert len(pc["ms1"]) == 3
    assert pc["ms1"][1][0] - pc["ms1"][0][0] == pytest.approx(1.00336, abs=1e-3)
    assert pc["ms1"][1][1] / pc["ms1"][0][1] == pytest.approx(0.45, abs=0.02)
    assert pc["ms2"] == [[float(mz), float(i)] for mz, i in zip(*spectra[pc["feature_id"]])]


# ---- SIRIUS (fake client) ------------------------------------------------------

def fake_models():
    """Stand-ins for the PySirius model classes: they keep their arguments."""
    def model(name):
        return lambda **kwargs: SimpleNamespace(model=name, **kwargs)
    names = ["BasicSpectrum", "SimplePeak", "FeatureImport", "JobSubmission", "Sirius", "Zodiac",
             "FingerprintPrediction", "Canopus", "StructureDbSearch", "MsNovelist", "SpectralLibrarySearch",
             "AccountCredentials"]
    return SimpleNamespace(**{name: model(name) for name in names})


class FakeSirius:
    """The part of the PySirius API used by atlas_ms.annotation.sirius, with canned results."""

    def __init__(self, logged_in=True, final_state="DONE", password="secret"):
        self.logged_in, self.final_state, self.password = logged_in, final_state, password
        self.imported, self.submission, self.polls, self.closed, self.logins = [], None, 0, [], []

    def account(self):
        def login(accept_terms, account_credentials):
            self.logins.append((accept_terms, account_credentials.username))
            if account_credentials.password != self.password:
                raise ValueError("(401) Unauthorized")
            self.logged_in = True
        return SimpleNamespace(is_logged_in=lambda: self.logged_in, login=login)

    def projects(self):
        return SimpleNamespace(get_projects=lambda: [], close_project=self.closed.append,
                               create_project=lambda pid, path_to_project: open(path_to_project, "w").close())

    def jobs(self):
        def start_job(pid, submission):
            self.submission = submission
            return SimpleNamespace(id="job-1")

        def get_job(pid, job_id):
            self.polls += 1  # running at the first poll, finished at the second
            state = "RUNNING" if self.polls == 1 else self.final_state
            return SimpleNamespace(progress=SimpleNamespace(state=state, current_progress=1, max_progress=2,
                                                            message="", error_message="boom"))
        return SimpleNamespace(start_job=start_job, get_job=get_job)

    def features(self):
        def formula(formula, zodiac, lipid=None, rank=1):
            return SimpleNamespace(formula_id=f"f{rank}", molecular_formula=formula, adduct="[M+H]+", zodiac_score=zodiac,
                                   sirius_score=100.0, sirius_score_normalized=0.99, num_of_explained_peaks=3,
                                   num_of_explainable_peaks=3, lipid_annotation=lipid)

        def get_formula_candidates(pid, afid, opt_fields):
            if afid == "af-1":  # a PC: El Gordo recognises the lipid
                lipid = SimpleNamespace(lipid_species="PC 34:1", lipid_class_name="Glycerophosphocholines",
                                        lipid_maps_id="LMGP01010005", hypothetical_structure="CCC")
                return [formula("C42H82NO8P", 0.95, lipid), formula("C40H80N4O7P", 0.03, rank=2)]
            return [formula("C10H12N2O", 0.40)]

        def get_structure_candidates(pid, afid, opt_fields):
            links = [SimpleNamespace(name="HMDB", id="1"), SimpleNamespace(name="LIPIDMAPS", id="2")]
            return [SimpleNamespace(structure_name="POPC", inchi_key="WTJKGGKOPKCXLL", smiles="CCC", molecular_formula="C42H82NO8P",
                                    adduct="[M+H]+", csi_score=-12.5, tanimoto_similarity=0.9, db_links=links)]

        def get_best_matching_compound_classes(pid, afid, formula_id):
            def cls(name, p):
                return SimpleNamespace(name=name, probability=p)
            return SimpleNamespace(classy_fire_lineage=[cls("Organic compounds", 1.0), cls("Lipids", 0.99),
                                                        cls("Phosphatidylcholines", 0.97)],
                                   npc_pathway=cls("Fatty acids", 0.9), npc_superclass=None, npc_class=cls("Glycerophosphocholines", 0.9))

        def get_de_novo_structure_candidates(pid, afid, opt_fields):
            return [SimpleNamespace(structure_name=None, inchi_key="NEWSTRUCTUREXX", smiles="CCCN", molecular_formula="C42H82NO8P",
                                    adduct="[M+H]+", csi_score=-20.0, tanimoto_similarity=0.7, db_links=[])]

        return SimpleNamespace(
            get_de_novo_structure_candidates=get_de_novo_structure_candidates,
            add_aligned_features=lambda pid, imports, profile: self.imported.extend(imports),
            get_aligned_features=lambda pid, opt_fields: [
                SimpleNamespace(external_feature_id=i.external_feature_id, aligned_feature_id=f"af-{i.external_feature_id}",
                                top_annotations=SimpleNamespace(confidence_exact_match=0.8, confidence_approx_match=0.9))
                for i in self.imported],
            get_formula_candidates=get_formula_candidates,
            get_structure_candidates=get_structure_candidates,
            get_best_matching_compound_classes=get_best_matching_compound_classes,
        )


ENTRIES = [
    {"feature_id": 1, "mz": 760.5851, "charge": 1, "rt_s": 180.0, "adducts": ["[M+H]+"],
     "ms1": [[760.5851, 1e6], [761.5885, 4.5e5]], "ms2": [[184.0733, 1e5], [104.107, 5e3]]},
    {"feature_id": 2, "mz": 177.1022, "charge": 1, "rt_s": 60.0, "adducts": [], "ms1": [], "ms2": [[120.08, 10.0]]},
]


def test_sirius_candidates_from_a_fake_sirius(tmp_path):
    api = FakeSirius()
    candidates = annotate_with_sirius(api, fake_models(), ENTRIES, SETTINGS, ["[M+H]+", "[M+Na]+"],
                                      tmp_path / "project.sirius", poll_s=0)
    # Input: our feature ids as external ids, isotope pattern only where measured.
    assert [i.external_feature_id for i in api.imported] == ["1", "2"]
    assert api.imported[0].merged_ms1.peaks[1].intensity == 4.5e5 and api.imported[1].merged_ms1 is None
    assert api.imported[0].detected_adducts == ["[M+H]+"] and api.imported[1].detected_adducts is None
    job = api.submission
    assert job.fallback_adducts == ["[M+H]+", "[M+Na]+"] and job.formula_id_params.profile == "ORBITRAP"
    assert job.ms_novelist_params.enabled and job.ms_novelist_params.number_of_candidate_to_predict == 64
    assert job.structure_db_search_params.structure_search_dbs == ["BIO"]
    assert api.polls == 2 and len(api.closed) == 1  # waited for the job; project closed

    rows = candidates.set_index(["feature_id", "source", "rank"])
    # Formulas: the top one is level 4 (ZODIAC 0.95), the second 5; feature 2's top one only 5 (ZODIAC 0.40).
    assert rows.loc[(1, "sirius:formula", 1), "proposed_level"] == "4"
    assert rows.loc[(1, "sirius:formula", 2), "proposed_level"] == "5"
    assert rows.loc[(2, "sirius:formula", 1), "proposed_level"] == "5"
    # El Gordo: the lipid species, level 3.
    assert rows.loc[(1, "sirius:elgordo", 1), ["name", "proposed_level"]].tolist() == ["PC 34:1", "3"]
    # CSI:FingerID: level 3, with the COSMIC confidence and the databases.
    csi = rows.loc[(1, "sirius:csi", 1)]
    assert (csi["name"], csi["proposed_level"]) == ("POPC", "3")
    assert json.loads(csi["evidence"])["cosmic_confidence_exact"] == 0.8
    assert json.loads(csi["evidence"])["databases"] == ["HMDB", "LIPIDMAPS"]
    # MSNovelist: de novo structures (named by their InChIKey), level 3, without COSMIC confidence.
    de_novo = rows.loc[(1, "sirius:msnovelist", 1)]
    assert (de_novo["name"], de_novo["proposed_level"]) == ("NEWSTRUCTUREXX", "3")
    assert "cosmic_confidence_exact" not in json.loads(de_novo["evidence"])
    # CANOPUS: the most specific ClassyFire class.
    canopus = rows.loc[(1, "sirius:canopus", 1)]
    assert canopus["name"] == "Phosphatidylcholines"
    assert json.loads(canopus["evidence"])["classyfire"] == "Organic compounds > Lipids > Phosphatidylcholines"


def test_sirius_login_and_failed_jobs(tmp_path):
    def run(api, account):
        return annotate_with_sirius(api, fake_models(), ENTRIES, SETTINGS, [], tmp_path / "p.sirius", poll_s=0,
                                    account=account)

    # Not logged in, no saved account: a clear message.
    with pytest.raises(RuntimeError, match="not logged in"):
        run(FakeSirius(logged_in=False), None)
    # The saved account needs the terms of service accepted.
    with pytest.raises(RuntimeError, match="terms of service"):
        run(FakeSirius(logged_in=False), {"username": "me@lab.org", "password": "secret", "accept_terms": False})
    # Wrong password: SIRIUS's message is passed on.
    with pytest.raises(RuntimeError, match="login failed for me@lab.org.*Unauthorized"):
        run(FakeSirius(logged_in=False), {"username": "me@lab.org", "password": "oops", "accept_terms": True})
    # Right password: logged in once, then the run goes on.
    api = FakeSirius(logged_in=False)
    assert not run(api, {"username": "me@lab.org", "password": "secret", "accept_terms": True}).empty
    assert api.logins == [(True, "me@lab.org")]
    # Already logged in (e.g. in the SIRIUS GUI): the saved account is not used.
    api = FakeSirius(logged_in=True)
    run(api, {"username": "me@lab.org", "password": "secret", "accept_terms": True})
    assert api.logins == []
    with pytest.raises(RuntimeError, match="FAILED: boom"):
        wait_for(FakeSirius(final_state="FAILED"), "pid", SimpleNamespace(id="job-1"), poll_s=0)


def test_sirius_account_file_is_private(tmp_path):
    import stat

    from atlas_ms.credentials import forget_sirius_account, load_sirius_account, save_sirius_account

    path = tmp_path / "config" / "sirius_account.json"
    save_sirius_account("me@lab.org", "secret", True, path)
    assert stat.S_IMODE(path.stat().st_mode) == 0o600  # readable by its owner only
    assert load_sirius_account(path) == {"username": "me@lab.org", "password": "secret", "accept_terms": True}
    save_sirius_account("me@lab.org", "new", True, path)  # overwriting keeps it private
    assert stat.S_IMODE(path.stat().st_mode) == 0o600 and load_sirius_account(path)["password"] == "new"
    forget_sirius_account(path)
    assert load_sirius_account(path) is None


def test_best_annotation_prefers_structures_to_classes(tmp_path):
    """Same level (3): CSI:FingerID's structure comes before CANOPUS's class (source priority)."""
    from atlas_ms.annotation.harmonize import harmonize
    from atlas_ms.config import HarmonizationSettings

    candidates = annotate_with_sirius(FakeSirius(), fake_models(), ENTRIES[1:], SETTINGS, [],
                                      tmp_path / "project.sirius", poll_s=0)
    features = pd.DataFrame({"feature_id": [2], "rt": [60.0]})
    nodes = pd.DataFrame({"feature_id": [2], "family": [-1], "qc": [""]})
    _, best = harmonize(candidates, features, nodes, HarmonizationSettings())
    assert best.iloc[0]["source"] == "sirius:csi"


# ---- MS2Query -----------------------------------------------------------------------

def test_ms2query_results(tmp_path):
    columns = ["query_spectrum_nr", "ms2query_model_prediction", "precursor_mz_difference", "precursor_mz_query_spectrum",
               "precursor_mz_analog", "inchikey", "analog_compound_name", "smiles", "scans", "cf_kingdom",
               "cf_superclass", "cf_class", "cf_subclass", "cf_direct_parent", "npc_class_results",
               "npc_superclass_results", "npc_pathway_results"]
    rows = [
        [1, 0.92, 0.0004, 760.5851, 760.5847, "WTJKGGKOPKCXLL", "PC(16:0/18:1)", "CCC", 3, "Organic compounds",
         "Lipids and lipid-like molecules", "Glycerophospholipids", "Glycerophosphocholines", "Phosphatidylcholines",
         "Glycerophosphocholines", "Glycerophospholipids", "Fatty acids"],
        [1, 0.75, 14.0157, 760.5851, 746.5694, "ABCDEFGHIJKLMN", "PC(15:0/18:1)", "CC", 3, None, None, None, None, None,
         None, None, None],
        [2, 0.71, 0.0, 496.3398, 496.3398, "ASWBNKHCZGQVJV", "LPC(16:0)", "C", 1, None, None, None, None, None, None,
         None, None],
    ]
    csv = tmp_path / "ms2_spectra.csv"
    pd.DataFrame(rows, columns=columns).to_csv(csv, index=False)
    candidates = parse_results(csv, precursor_tolerance_ppm=10.0).set_index(["feature_id", "rank"])
    best = candidates.loc[(3, 1)]
    assert (best["name"], best["source"], best["proposed_level"]) == ("PC(16:0/18:1)", "ms2query", "3")
    evidence = json.loads(best["evidence"])
    assert evidence["exact_precursor_match"] is True and evidence["classyfire"].endswith("Phosphatidylcholines")
    assert json.loads(candidates.loc[(3, 2), "evidence"])["exact_precursor_match"] is False
    assert candidates.loc[(1, 1), "name"] == "LPC(16:0)"


# ---- Workflow ---------------------------------------------------------------------------

def test_sirius_and_ms2query_join_the_workflow_when_switched_on(processed_project, tmp_path):
    root = tmp_path / "copy"
    shutil.copytree(processed_project.root, root, symlinks=True)
    config = yaml.safe_load((root / "project.yaml").read_text())
    config["sirius"]["enabled"] = True
    config["ms2query"]["enabled"] = True
    config["ms2query"]["models_dir"] = str(tmp_path)  # a models folder: no download step
    (root / "project.yaml").write_text(yaml.safe_dump(config))
    jobs = planned_jobs(snakemake(root, "--dry-run"))
    assert jobs == {"prepare_sirius_input": 1, "run_sirius": 1, "run_ms2query": 1, "harmonize": 1,
                    "export_graphml": 1, "all": 1}
