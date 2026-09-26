"""The app, without a browser: views built from the processed synthetic study."""

import asyncio
import shutil
import socket
import sys

import numpy as np
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
    assert [name for name in app.tabs._names] == ["Setup", "Network", "Annotation"]
    network = app.network
    features = network.table.value
    network.selected = [5]
    assert features.iloc[network.table.selection]["feature_id"].tolist() == [5]
    network.table.selection = [0, 2]
    assert network.selected == features.iloc[[0, 2]]["feature_id"].tolist()
    # One selected feature: header, MS2 spectrum, chromatograms.
    network.selected = [5]
    assert len(network.details().objects) == 3
    # Two selected features with MS2: header, mirror caption, mirror plot, chromatograms.
    network.selected = [1, 3]
    assert len(network.details().objects) == 4
    # No MS1 scan around the feature's RT (empty chromatograms): a warning instead of the plot.
    assert network.reader.xic(features["mz"].iloc[0], 1e5).empty
    network.reader.xic = lambda mz, rt: pd.DataFrame(columns=["sample", "rt", "intensity"])
    network.selected = [5]
    assert type(network.details().objects[-1]).__name__ == "Alert"
    # Every colouring and every node size (one per sample) can be computed,
    # and with every size, no two nodes overlap.
    coords = network.nodes[["x", "y"]].to_numpy()
    dist = np.sqrt(((coords[:, None] - coords[None]) ** 2).sum(-1))
    np.fill_diagonal(dist, np.inf)
    for color in network.param.color_by.objects:
        for size in network.param.size_by.objects:
            network._points(color, size)
            radii = network._radii(size)
            assert (dist > radii[:, None] + radii[None, :]).all()
    assert "treat_2" in network.param.size_by.objects


def test_feature_table_in_minutes_with_search_boxes(processed_project):
    network = AtlasApp(str(processed_project.root)).network
    table = network.table
    assert "rt_min" in table.value and "rt" not in table.value
    cer = table.value[table.value["mz"].round(2) == 538.52]
    assert cer["rt_min"].iloc[0] == pytest.approx(200 / 60, abs=0.01)
    assert set(table.header_filters) == set(table.value.columns)
    assert table.page_size == 25


def test_fragment_and_neutral_loss_search(processed_project):
    network = AtlasApp(str(processed_project.root)).network
    features = network.results.features
    ids = {name: feature_of(features, name)["feature_id"] for name in (
        "PC 34:1 [M+H]+", "LPC 16:0 [M+H]+", "Cer 34:1;O2 [M+H]+")}

    # Phosphocholine fragments: PC and LPC have both, SM only m/z 184.
    network.fragment_text.value = "184.0733, 104.107"
    network.search_fragments()
    assert sorted(network.matched) == sorted([ids["PC 34:1 [M+H]+"], ids["LPC 16:0 [M+H]+"]])

    # Water loss: LPC (496.34 -> 478.33) and Cer (538.52 -> 520.51).
    network.fragment_kind.value = "neutral losses"
    network.fragment_text.value = "18.0106"
    network.search_fragments()
    assert sorted(network.matched) == sorted([ids["LPC 16:0 [M+H]+"], ids["Cer 34:1;O2 [M+H]+"]])

    network.clear_search()
    assert network.matched == []


def test_peak_labels_do_not_overlap():
    import numpy as np

    from atlas_ms.app.plots import label_peaks

    mz = np.array([184.07, 184.5, 300.0, 478.3, 104.1])
    intensity = np.array([100.0, 90.0, 20.0, 10.0, 1.0])
    chosen = label_peaks(mz, intensity)
    # 184.5 is too close to the more intense 184.07; 104.1 is below 5% of the base peak.
    assert sorted(mz[chosen]) == [184.07, 300.0, 478.3]


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


def test_annotation_tab_follows_the_selection(processed_project):
    from atlas_ms.app.annotation_view import evidence_marks

    app = AtlasApp(str(processed_project.root))
    network, annotation = app.network, app.annotation
    features = network.results.features
    pc = int(feature_of(features, "PC 34:1 [M+H]+")["feature_id"])

    # The network table shows the best annotation and its level.
    row = network.table.value.set_index("feature_id").loc[pc]
    assert (row["annotation"], row["level"], row["lipid_class"]) == ("PC 16:0_18:1", "1", "PC")

    network.selected = [pc]
    candidates = annotation.table.value
    # Most confident first: the reference standard (L1), then the lipid rule (L3).
    assert candidates[["level", "source"]].values.tolist()[:2] == [["1", "library:standards"], ["3", "lipid_rules"]]
    assert annotation.table.selection == [0]
    # Library hit: a mirror plot of the feature against the library spectrum.
    plot = annotation.evidence(annotation.table.value, [0])
    assert type(plot).__name__ == "Bokeh"
    # Lipid rule: the spectrum with the diagnostic ions marked.
    lipid_rule = annotation._rows.iloc[1]
    assert [text for _, text in evidence_marks(lipid_rule["evidence"], 760.5851)] == ["184.0733", "104.1070"]
    assert type(annotation.evidence(annotation.table.value, [1])).__name__ == "Bokeh"
    assert "L1 · PC 16:0_18:1" in annotation.header(network.selected).object

    # The family lists PC's family members; clicking one selects it everywhere.
    members = annotation.family_table.value
    lpc = int(feature_of(features, "LPC 16:0 [M+H]+")["feature_id"])
    assert pc in set(members["feature_id"]) and lpc in set(members["feature_id"])
    annotation.family_table.selection = [int(members.index[members["feature_id"] == lpc][0])]
    assert network.selected == [lpc]
    assert annotation.table.value["name"].iloc[0] == "LPC 16:0"


def test_level_and_class_colours_have_a_legend(processed_project):
    network = AtlasApp(str(processed_project.root)).network
    network.color_by = "confidence level"
    assert "L1" in network.legend().object and "not annotated" in network.legend().object
    network.color_by = "lipid class"
    assert "TG" in network.legend().object
    network.color_by = "family"
    assert network.legend().object == ""


def test_setup_edits_the_library_list(project_copy, tmp_path):
    setup = SetupTab(project_copy)
    assert list(setup.libraries.value["name"]) == ["standards", "insilico"]
    setup._add_library(None)
    assert not setup.save()  # the new library has no file yet
    library = tmp_path / "extra.mgf"
    library.write_text("BEGIN IONS\nPEPMASS=500.0\n100.0 10.0\nEND IONS\n")
    setup.libraries.value.loc[2, ["name", "path"]] = ["extra", str(library)]
    assert setup.save()
    saved = project_copy.load_config().library_search
    assert [lib["name"] for lib in saved.libraries] == ["standards", "insilico", "extra"]
    assert saved.library("extra")["path"] == str(library.resolve())
