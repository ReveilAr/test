"""The app, without a browser: views built from the processed synthetic study."""

import asyncio
import json
import shutil
import socket
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from atlas_ms.app.data import load_results, read_mgf
from atlas_ms.app.main import AtlasApp, free_port
from atlas_ms.app.run_view import RunPanel
from atlas_ms.app.setup_view import SetupTab
from atlas_ms.preprocessing.msdata import RunReader
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
    reader = RunReader(processed_project.root, list(processed_project.load_samples().index))
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


# ---- Drag and drop, SIRIUS account, Quit ----------------------------------------

def drop(dropper, name: str, content: bytes, pieces: int = 2) -> None:
    """Simulate the browser uploading a file to a drop area, in ``pieces`` pieces."""
    from types import SimpleNamespace

    size = -(-len(content) // pieces)
    for i in range(pieces):
        dropper._process_event(SimpleNamespace(event_name="upload_event", data={
            "name": name, "chunk": i + 1, "total_chunks": pieces, "type": "", "data": content[i * size:(i + 1) * size]}))


def test_dropped_files_are_written_piece_by_piece(tmp_path):
    from atlas_ms.app.uploads import drop_area

    dropper = drop_area(str(tmp_path / "raw"), [".raw", ".mzml"], "drop")
    drop(dropper, "run_1.RAW", b"0123456789", pieces=3)
    assert (tmp_path / "raw" / "run_1.RAW").read_bytes() == b"0123456789"
    assert dropper.saved == [str(tmp_path / "raw" / "run_1.RAW")] and not list((tmp_path / "raw").glob("*.part"))
    dropper.get_root()  # renders in a page (only FileDropper's own settings reach the browser)
    drop(dropper, "notes.txt", b"x")  # wrong type
    drop(dropper, "../../escape.mzML", b"x")  # only the name is kept
    assert dropper.rejected == ["notes.txt"] and (tmp_path / "raw" / "escape.mzML").exists()


def test_setup_drops_samples_libraries_and_models(project_copy, tmp_path):
    from atlas_ms.app.setup_view import MS2QUERY_DROP_DIR

    setup = SetupTab(project_copy, account_file=tmp_path / "account.json")
    drop(setup.sample_drop, "new_run.mzML", b"<mzML/>")
    assert setup.samples.value["sample"].tolist()[-1] == "new_run"
    drop(setup.library_drop, "My lib.msp", b"NAME: x\n")
    library = setup.libraries.value.iloc[-1]
    assert (library["name"], library["kind"]) == ("My_lib", "experimental")
    assert Path(library["path"]) == project_copy.root / "libraries" / "My lib.msp"
    setup.ms2deepscore_drop.folder = str(tmp_path / "models")  # not the real model cache
    drop(setup.ms2deepscore_drop, "custom.pt", b"weights")
    assert setup.config.scoring.ms2deepscore_model == str(tmp_path / "models" / "custom.pt")
    setup.ms2query_drop.folder = str(tmp_path / "ms2query")
    drop(setup.ms2query_drop, "library.sqlite", b"db")
    assert setup.config.ms2query.models_dir == str(MS2QUERY_DROP_DIR)
    # Remove the dropped sample again, and save.
    setup.samples.selection = [len(setup.samples.value) - 1]
    setup._remove_samples(None)
    assert setup.save()
    assert "new_run" not in project_copy.load_samples().index


def test_sirius_login_is_saved_outside_the_project(project_copy, tmp_path):
    account_file = tmp_path / "account.json"
    setup = SetupTab(project_copy, account_file=account_file)
    setup.sirius_user.value, setup.sirius_password.value = "me@lab.org", "secret"
    assert not setup.save_sirius_login()  # terms of service not accepted
    setup.sirius_terms.value = True
    assert setup.save_sirius_login()
    assert json.loads(account_file.read_text())["username"] == "me@lab.org"
    assert setup.sirius_password.value == "" and "me@lab.org" in setup.sirius_status.object
    assert setup.save() and "secret" not in project_copy.config_path.read_text()
    setup.forget_sirius_login()
    assert not account_file.exists()


def test_new_project_from_dropped_files(tmp_path, study_files):
    app = AtlasApp()
    app.path.value = str(tmp_path / "dropped_project")
    for file in study_files[:2]:
        drop(app.raw_drop, file.name, file.read_bytes(), pieces=4)
    app.create()
    samples = Project(tmp_path / "dropped_project").load_samples()
    assert list(samples.index) == ["ctrl_1", "ctrl_2"]
    assert Path(samples.loc["ctrl_1", "file"]).parent == tmp_path / "dropped_project" / "raw"


def test_quit_asks_twice_then_stops(processed_project):
    stopped = []
    app = AtlasApp(str(processed_project.root), on_quit=lambda: stopped.append(True))
    assert app.quit_button.visible
    app.quit()
    assert not stopped and app.quit_button.label == "Click again to quit"
    app.quit()
    assert stopped == [True] and "stopped" in app.content.objects[0].object
    assert not AtlasApp(str(processed_project.root)).quit_button.visible  # no server to stop


def test_stopping_the_server_frees_its_port():
    """A real server, stopped as by the Quit button: the process ends and the port is free."""
    port = free_port(5390)
    script = ("import threading; import panel as pn; from atlas_ms.app.main import stop_later\n"
              f"server = pn.serve(lambda: pn.pane.Markdown('hi'), port={port}, show=False, start=False)\n"
              "server.start(); stop_later(server, 0.5); server.io_loop.start(); print('stopped')")
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0 and "stopped" in result.stdout, result.stderr[-2000:]
    assert free_port(port) == port


def test_runner_limits_memory_and_threads(tmp_path):
    from atlas_ms.runner import available_cores, snakemake_command, snakemake_env

    command = snakemake_command(tmp_path)
    assert command[command.index("--cores") + 1] == str(available_cores())
    memory = int(command[command.index("--resources") + 1].split("=")[1])
    assert 0 < memory < 10_000_000
    assert snakemake_env()["OMP_NUM_THREADS"] == "1"
