"""
End-to-end tests of the preprocessing workflow on the synthetic study
(see synthetic.py): real Snakemake runs, real pyOpenMS processing.
"""

import re
import shutil
from pathlib import Path

import pandas as pd
import pyopenms as oms
import pytest
import yaml

import atlas_ms
from atlas_ms.config import ExportSettings
from atlas_ms.preprocessing.export import export_gnps
from atlas_ms.preprocessing.msdata import load_consensus_map
from atlas_ms.project import Project
from conftest import snakemake
from synthetic import LIPIDS


def planned_jobs(result) -> dict[str, int]:
    """Rule name -> job count, from the "Job stats" table of a Snakemake (dry) run."""
    text = result.stdout + result.stderr
    if "Job stats:" not in text:
        return {}
    table = text.split("Job stats:", 1)[1].split("total", 1)[0]
    return {m.group(1): int(m.group(2)) for m in re.finditer(r"^(\w+)\s+(\d+)\s*$", table, re.MULTILINE)}


def feature_of(features: pd.DataFrame, name: str) -> pd.Series:
    """The feature matching a synthetic compound (m/z within 5 ppm, RT within 5 s)."""
    compound = next(c for c in LIPIDS if c.name == name)
    match = features[
        ((features["mz"] - compound.mz).abs() / compound.mz < 5e-6)
        & ((features["rt"] - compound.rt).abs() < 5)
    ]
    assert len(match) == 1, f"{name}: {len(match)} matching features"
    return match.iloc[0]


def read_mgf_headers(path: Path) -> pd.DataFrame:
    """SCANS and PEPMASS of every spectrum of an MGF file."""
    blocks = path.read_text().split("BEGIN IONS")[1:]
    rows = []
    for block in blocks:
        fields = dict(line.split("=", 1) for line in block.splitlines() if "=" in line)
        rows.append({"scans": int(fields["SCANS"]), "pepmass": float(fields["PEPMASS"])})
    return pd.DataFrame(rows)


@pytest.fixture(scope="module")
def results(processed_project):
    root = processed_project.results_dir
    return {
        "features": pd.read_parquet(root / "features.parquet"),
        "quant": pd.read_parquet(root / "quant.parquet").set_index("feature_id"),
        "gnps": root / "gnps",
    }


def test_every_compound_becomes_one_feature(results):
    features = results["features"]
    assert len(features) == len(LIPIDS)
    for compound in LIPIDS:
        feature_of(features, compound.name)


def test_features_with_ms2_come_first(results):
    features = results["features"]
    assert list(features["feature_id"]) == list(range(1, len(features) + 1))
    # Once a feature without MS2 appears, no feature with MS2 may follow.
    has_ms2 = list(features["has_ms2"])
    assert has_ms2 == sorted(has_ms2, reverse=True)
    assert not feature_of(features, "DG 34:1 [M+NH4]+")["has_ms2"]


def test_gap_filling_recovers_the_weak_compound(results):
    cer = feature_of(results["features"], "Cer 34:1;O2 [M+H]+")
    values = results["quant"].loc[cer["feature_id"]]
    assert values.notna().all()
    # In treat_2 the compound was simulated 300 times weaker.
    ratio = values["ctrl_1"] / values["treat_2"]
    assert 150 < ratio < 600


def test_quantification_follows_the_simulated_change(results):
    tg = feature_of(results["features"], "TG 52:2 [M+NH4]+")
    values = results["quant"].loc[tg["feature_id"]]
    fold_change = values[["treat_1", "treat_2"]].mean() / values[["ctrl_1", "ctrl_2"]].mean()
    assert 2.5 < fold_change < 3.5


def test_adducts_of_one_molecule_are_grouped(results):
    features = results["features"]
    protonated = feature_of(features, "PC 34:1 [M+H]+")
    sodiated = feature_of(features, "PC 34:1 [M+Na]+")
    assert protonated["ion"] == "[M+H]+"
    assert sodiated["ion"] == "[M+Na]+"
    assert protonated["ion_network"] == sodiated["ion_network"]
    assert protonated["ion_partners"] == str(sodiated["feature_id"])
    # Features without a partner get no adduct: it would only be a guess.
    assert features["ion"].notna().sum() == 2


def test_gnps_export_uses_the_same_feature_ids(results):
    features = results["features"].set_index("feature_id")
    gnps = results["gnps"]

    mgf = read_mgf_headers(gnps / "ms2_spectra.mgf")
    ms2_ids = set(features.index[features["has_ms2"]])
    assert set(mgf["scans"]) == ms2_ids
    for row in mgf.itertuples():
        assert row.pepmass == pytest.approx(features.loc[row.scans, "mz"], rel=1e-6)

    quant_rows = [line.split("\t") for line in (gnps / "quantification_table.txt").read_text().splitlines()]
    header = next(r for r in quant_rows if r[0] == "#CONSENSUS")
    row_ids = {int(r[header.index("row ID")]) for r in quant_rows if r[0] == "CONSENSUS"}
    assert row_ids == ms2_ids

    # The metadata must name the same mzML files as the quantification table.
    map_files = {Path(r[2]).name for r in quant_rows if r[0] == "MAP"}
    metadata = pd.read_csv(gnps / "metadata.tsv", sep="\t")
    assert set(metadata["filename"]) == map_files
    assert list(metadata["ATTRIBUTE_group"]) == ["ctrl", "ctrl", "treat", "treat"]

    pairs = pd.read_csv(gnps / "iimn_supplementary_pairs.csv")
    pc_ids = {feature_of(results["features"], n)["feature_id"] for n in ("PC 34:1 [M+H]+", "PC 34:1 [M+Na]+")}
    assert {int(pairs.iloc[0]["ID1"]), int(pairs.iloc[0]["ID2"])} == pc_ids


def test_only_affected_steps_are_rerun(processed_project, tmp_path):
    # Work on a copy: the processed project is shared by other tests.
    # copytree keeps file times (and the mzML symbolic links).
    root = tmp_path / "copy"
    shutil.copytree(processed_project.root, root, symlinks=True)
    assert planned_jobs(snakemake(root, "--dry-run")) == {}  # nothing to do

    config = yaml.safe_load((root / "project.yaml").read_text())
    config["export"]["min_detection_fraction"] = 0.5
    (root / "project.yaml").write_text(yaml.safe_dump(config))
    assert planned_jobs(snakemake(root, "--dry-run")) == {"export": 1, "all": 1}

    config["linking"]["rt_tol_s"] = 20.0
    (root / "project.yaml").write_text(yaml.safe_dump(config))
    jobs = planned_jobs(snakemake(root, "--dry-run"))
    assert "find_features" not in jobs and "align" not in jobs
    assert jobs["link"] == 1 and jobs["export"] == 1


def test_thermo_files_are_converted_and_gap_filling_can_be_skipped(tmp_path):
    raw = tmp_path / "run1.raw"
    raw.touch()
    project = Project.create(tmp_path / "proj", [raw])
    config = project.load_config()
    config.gap_filling.enabled = False
    project.save_config(config)
    jobs = planned_jobs(snakemake(project.root, "--dry-run"))
    assert jobs["convert_thermo"] == 1
    assert "link_mzml" not in jobs
    assert "fill_gaps" not in jobs and "plan_gap_filling" not in jobs


def test_gnps_export_without_any_ms2(processed_project, tmp_path):
    """An MS1-only study must still produce every GNPS file (empty ones)."""
    root = processed_project.root
    consensus = load_consensus_map(root / "work/consensus/gap_filled.consensusXML")
    no_ms2 = oms.ConsensusMap(consensus)
    no_ms2.clear(False)
    for cf in consensus:
        cf.setPeptideIdentifications(oms.PeptideIdentificationList())
        no_ms2.push_back(cf)
    samples = processed_project.load_samples()
    mzml = [root / "work/mzml" / f"{name}.mzML" for name in samples.index]
    export_gnps(no_ms2, mzml, samples, tmp_path / "gnps", tmp_path / "gnps.consensusXML", ExportSettings())
    for name in ("ms2_spectra.mgf", "quantification_table.txt", "metadata.tsv", "iimn_supplementary_pairs.csv"):
        assert (tmp_path / "gnps" / name).exists(), name


def test_rule_logs_record_the_code_version(processed_project):
    """Every rule log starts with the version and commit that produced it."""
    first_line = (processed_project.root / "logs" / "export.log").read_text().splitlines()[0]
    assert f"ATLAS-MS {atlas_ms.__version__}, commit " in first_line
