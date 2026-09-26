"""Network construction steps (pure logic) and the network of the synthetic study."""

import shutil

import networkx as nx
import numpy as np
import pandas as pd
import pytest
import yaml
from ms2deepscore.models import SiameseSpectralModel
from ms2deepscore.SettingsMS2Deepscore import SettingsMS2Deepscore

from atlas_ms.config import NetworkSettings, ProjectConfig, ScoringSettings
from atlas_ms.network.graph import filter_edges, limit_family_size, number_groups
from atlas_ms.network.model_files import DEFAULT_MODEL_PATH
from atlas_ms.network.scoring import candidate_edges, score_spectra
from atlas_ms.network.spectra import load_spectra
from conftest import snakemake
from test_workflow import feature_of, planned_jobs


def tiny_model(path) -> str:
    """An untrained, very small MS2DeepScore model (the real one is downloaded on first use)."""
    SiameseSpectralModel(SettingsMS2Deepscore(base_dims=(32,), embedding_dim=8)).save(str(path))
    return str(path)


def edges(rows):
    return pd.DataFrame(rows, columns=["source", "target", "score", "matched_peaks"])


def test_candidates_keep_the_best_neighbours_once():
    scores = np.array([
        [1.0, 0.9, 0.8, 0.1],
        [0.9, 1.0, 0.5, 0.2],
        [0.8, 0.5, 1.0, 0.4],
        [0.1, 0.2, 0.4, 1.0],
    ])
    settings = ScoringSettings(candidates_per_spectrum=1, min_candidate_score=0.3)
    pool = candidate_edges(scores, None, [10, 20, 30, 40], settings)
    # Best neighbour of 10 -> 20, of 20 -> 10 (same edge), of 30 -> 10, of 40 -> 30.
    assert sorted(zip(pool["source"], pool["target"])) == [(10, 20), (10, 30), (30, 40)]


def test_top_k_needs_both_ends():
    # Node 1 is a hub: every other node's best neighbour, but with K=2 it
    # keeps only its two best edges.
    candidates = edges([(1, 2, 0.95, 5), (1, 3, 0.90, 5), (1, 4, 0.85, 5), (3, 4, 0.80, 5)])
    kept = filter_edges(candidates, NetworkSettings(min_score=0.7, min_matched_peaks=1, top_k=2))
    assert sorted(zip(kept["source"], kept["target"])) == [(1, 2), (1, 3), (3, 4)]


def test_cutoffs():
    candidates = edges([(1, 2, 0.95, 2), (1, 3, 0.65, 9), (2, 3, 0.90, 9)])
    kept = filter_edges(candidates, NetworkSettings(min_score=0.7, min_matched_peaks=3))
    assert list(zip(kept["source"], kept["target"])) == [(2, 3)]


def test_large_families_lose_their_weakest_edges():
    # A chain 1-2-3-4-5 whose middle edge is the weakest.
    chain = edges([(1, 2, 0.9, 5), (2, 3, 0.8, 5), (3, 4, 0.7, 5), (4, 5, 0.95, 5)])
    kept = limit_family_size(chain, max_size=3)
    assert sorted(zip(kept["source"], kept["target"])) == [(1, 2), (2, 3), (4, 5)]
    assert len(limit_family_size(chain, max_size=0)) == 4  # 0 = no limit


def test_groups_are_numbered_by_size():
    numbering = number_groups([{5}, {1, 2}, {3, 4, 6}])
    assert numbering == {3: 1, 4: 1, 6: 1, 1: 2, 2: 2}  # singleton 5 not numbered


def test_network_of_the_synthetic_study(processed_project):
    features = pd.read_parquet(processed_project.results_dir / "features.parquet")
    nodes = pd.read_parquet(processed_project.results_dir / "network" / "nodes.parquet").set_index("feature_id")
    network_edges = pd.read_parquet(processed_project.results_dir / "network" / "edges.parquet")

    # Nodes are the features with MS2: all lipids except the DG.
    assert set(nodes.index) == set(features.loc[features["has_ms2"], "feature_id"])

    ids = {name: feature_of(features, name)["feature_id"] for name in (
        "PC 34:1 [M+H]+", "PC 34:1 [M+Na]+", "SM 34:1;O2 [M+H]+", "LPC 16:0 [M+H]+",
        "Cer 34:1;O2 [M+H]+", "TG 52:2 [M+NH4]+", "CE 18:1 [M+NH4]+")}
    # The choline lipids share the m/z 184 fragment: one family, with the
    # sodiated PC joined by its ion-identity edge.
    choline = [ids[n] for n in ("PC 34:1 [M+H]+", "SM 34:1;O2 [M+H]+", "LPC 16:0 [M+H]+", "PC 34:1 [M+Na]+")]
    assert nodes.loc[choline, "family"].nunique() == 1
    assert nodes.loc[choline[0], "family"] == 1
    for name in ("Cer 34:1;O2 [M+H]+", "TG 52:2 [M+NH4]+", "CE 18:1 [M+NH4]+"):
        assert nodes.loc[ids[name], "family"] == -1  # singletons

    ion = network_edges[network_edges["edge_type"] == "ion_identity"]
    assert {frozenset(p) for p in zip(ion["source"], ion["target"])} == {
        frozenset((ids["PC 34:1 [M+H]+"], ids["PC 34:1 [M+Na]+"]))}
    # PC [M+H]+ - SM and PC [M+H]+ - LPC share 2 fragments (directly or
    # shifted by the precursor difference); SM - LPC share only m/z 184.
    spectral = network_edges[network_edges["edge_type"] == "spectral"]
    assert {frozenset(p) for p in zip(spectral["source"], spectral["target"])} == {
        frozenset((ids["PC 34:1 [M+H]+"], ids["SM 34:1;O2 [M+H]+"])),
        frozenset((ids["PC 34:1 [M+H]+"], ids["LPC 16:0 [M+H]+"])),
    }
    assert (spectral["score"] >= 0.7).all() and (spectral["matched_peaks"] >= 2).all()

    # GraphML for Cytoscape, with the feature table as node attributes.
    graph = nx.read_graphml(processed_project.results_dir / "network" / "network.graphml")
    assert graph.number_of_nodes() == len(nodes) and graph.number_of_edges() == len(network_edges)
    assert "mz" in next(iter(graph.nodes(data=True)))[1]


def test_ms2deepscore_scoring(processed_project, tmp_path):
    """MS2DeepScore scoring, with a tiny untrained model."""
    spectra = load_spectra(processed_project.root / "work" / "network" / "spectra.pickle")
    settings = ScoringSettings(score="ms2deepscore", ms2deepscore_model=tiny_model(tmp_path / "tiny.pt"))
    pool = score_spectra(spectra, settings)
    # An untrained model scores everything alike: every pair is a candidate.
    n = len(spectra)
    assert len(pool) == n * (n - 1) // 2
    assert pool["matched_peaks"].isna().all()
    assert pool["score"].between(-1, 1).all()


def test_config_rejects_inconsistent_network_cutoffs():
    with pytest.raises(ValueError, match="min_candidate_score"):
        ProjectConfig.from_dict({"network": {"min_score": 0.2}})
    with pytest.raises(ValueError, match="candidates_per_spectrum"):
        ProjectConfig.from_dict({"network": {"top_k": 80}})


def test_ms2deepscore_through_the_workflow(processed_project, tmp_path):
    """Switching to MS2DeepScore re-runs only scoring and network; the model is a rule input."""
    root = tmp_path / "copy"
    shutil.copytree(processed_project.root, root, symlinks=True)
    config = yaml.safe_load((root / "project.yaml").read_text())

    # Without a model file, the pretrained model would be downloaded first.
    config["scoring"]["score"] = "ms2deepscore"
    (root / "project.yaml").write_text(yaml.safe_dump(config))
    jobs = planned_jobs(snakemake(root, "--dry-run"))
    if not DEFAULT_MODEL_PATH.exists():
        assert jobs["download_ms2deepscore_model"] == 1

    # With a local model: scoring and network only, and the run succeeds.
    config["scoring"]["ms2deepscore_model"] = tiny_model(tmp_path / "tiny.pt")
    config["network"]["min_score"] = 0.3  # an untrained model's scores are arbitrary
    (root / "project.yaml").write_text(yaml.safe_dump(config))
    assert planned_jobs(snakemake(root, "--dry-run")) == {"score_spectra": 1, "build_network": 1, "all": 1}
    result = snakemake(root)
    assert result.returncode == 0, result.stderr[-3000:]
    assert (root / "results" / "network" / "nodes.parquet").exists()



def test_layout_keeps_nodes_and_families_apart():
    from atlas_ms.network.graph import MIN_DISTANCE, layout

    # Three dense families (cliques of 30, 12 and 5 nodes) and some singletons.
    graph = nx.Graph()
    start = 0
    for size in (30, 12, 5):
        graph.add_edges_from((start + i, start + j) for i in range(size) for j in range(i + 1, size))
        start += size
    graph.add_nodes_from(range(start, start + 8))
    positions = layout(graph)
    assert set(positions) == set(graph)

    for family in nx.connected_components(graph):
        coords = np.array([positions[n] for n in family])
        if len(coords) > 1:
            dist = np.sqrt(((coords[:, None] - coords[None]) ** 2).sum(-1))
            np.fill_diagonal(dist, np.inf)
            assert dist.min() >= MIN_DISTANCE * 0.99  # no overlapping nodes inside a family
    # Nor anywhere else (between families, or with the singletons).
    coords = np.array(list(positions.values()))
    dist = np.sqrt(((coords[:, None] - coords[None]) ** 2).sum(-1))
    np.fill_diagonal(dist, np.inf)
    assert dist.min() >= MIN_DISTANCE * 0.99
    # The bounding boxes of different families do not overlap.
    boxes = [np.array([positions[n] for n in f]) for f in nx.connected_components(graph) if len(f) > 1]
    for i, a in enumerate(boxes):
        for b in boxes[i + 1:]:
            separate_x = a[:, 0].max() < b[:, 0].min() or b[:, 0].max() < a[:, 0].min()
            separate_y = a[:, 1].max() < b[:, 1].min() or b[:, 1].max() < a[:, 1].min()
            assert separate_x or separate_y
