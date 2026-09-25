"""The app, without a browser: views built from the processed synthetic study."""

import asyncio
import shutil
import socket
import sys

import pandas as pd
import pytest

from atlas_ms.app.data import ChromatogramReader, load_results, read_mgf
from atlas_ms.app.main import AtlasApp, free_port
from atlas_ms.app.run_view import RunPanel
from atlas_ms.app.setup_view import SetupTab
from atlas_ms.project import Project
from test_workflow import feature_of


@pytest.fixture
def project_copy(processed_project, tmp_path):
    """A copy of the processed project, for tests that save changes."""
    root = tmp_path / "copy"
    shutil.copytree(processed_project.root, root, symlinks=True)
    return Project(root)


def test_mgf_reader(processed_project):
    spectra = read_mgf(processed_project.results_dir / "gnps" / "ms2_spectra.mgf")
    features = pd.read_parquet(processed_project.results_dir / "features.parquet")
    assert set(spectra) == set(features.loc[features["has_ms2"], "feature_id"])
    lpc = feature_of(features, "LPC 16:0 [M+H]+")["feature_id"]
    assert 184.0733 in spectra[lpc][0].round(4)


def test_chromatograms_follow_the_simulated_intensities(processed_project):
    results = load_results(processed_project)
    cer = feature_of(results.features, "Cer 34:1;O2 [M+H]+")
    reader = ChromatogramReader(processed_project, list(processed_project.load_samples().index))
    apex = reader.xic(cer["mz"], cer["rt"]).groupby("sample")["intensity"].max()
    assert apex["ctrl_1"] / apex["treat_2"] == pytest.approx(300, rel=0.05)


def test_network_and_table_share_the_selection(processed_project):
    app = AtlasApp(str(processed_project.root))
    assert [name for name in app.tabs._names] == ["Setup", "Network"]
    network = app.network
    features = network.table.value
    network.selected = [5]
    assert features.iloc[network.table.selection]["feature_id"].tolist() == [5]
    network.table.selection = [0, 2]
    assert network.selected == features.iloc[[0, 2]]["feature_id"].tolist()
    # Details of the first selected feature: header, MS2 spectrum, chromatograms.
    assert len(network.details().objects) == 3
    for option in network.param.color_by.objects:
        network.color_by = option  # every colouring can be computed


def test_setup_saves_metadata_and_rejects_bad_adducts(project_copy):
    setup = SetupTab(project_copy)
    setup.new_column.value = "species"
    setup._add_column(None)
    setup.samples.value = setup.samples.value.assign(ATTRIBUTE_species=["a", "a", "b", "b"])
    setup.config.apply_instrument_preset("qtof")
    assert setup.save()
    assert list(project_copy.load_samples()["ATTRIBUTE_species"]) == ["a", "a", "b", "b"]
    assert project_copy.load_config().instrument.mass_error_ppm == 20.0

    # A negative-mode adduct in positive mode is refused, and nothing is written.
    setup.adducts.value = pd.DataFrame([{"name": "[M-H]-", "openms": "H-1:-", "probability": 1.0}])
    assert not setup.save()
    assert "Not saved" in setup.message.object
    assert project_copy.load_config().adducts.adducts[0]["name"] == "[M+H]+"


def test_run_panel_shows_progress_and_reloads(project_copy):
    fake_snakemake = "print('1 of 2 steps (50%) done'); print('2 of 2 steps (100%) done')"
    reloaded = []
    panel = RunPanel(
        project_copy, before_run=lambda: True, after_run=lambda: reloaded.append(True),
        command=lambda root, cores: [sys.executable, "-c", fake_snakemake],
    )
    assert asyncio.run(panel.run()) == 0
    assert panel.progress.value == 100 and reloaded == [True]
    assert "2 of 2 steps" in panel.log.object


def test_a_busy_port_is_skipped():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as busy:
        busy.bind(("", 0))  # any free port, now taken
        busy.listen()
        port = busy.getsockname()[1]
        assert free_port(port) != port
