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
from test_workflow import AFTER_NETWORK, feature_of, planned_jobs


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

    # GraphML for Cytoscape, with the feature table and the best annotation as node attributes.
    graph = nx.read_graphml(processed_project.results_dir / "network" / "network.graphml")
    assert graph.number_of_nodes() == len(nodes) and graph.number_of_edges() == len(network_edges)
    pc = graph.nodes[str(ids["PC 34:1 [M+H]+"])]
    assert pc["mz"] == pytest.approx(760.5851, abs=0.001)
    assert pc["annotation"] == "PC 16:0_18:1" and pc["level"] == "1" and pc["lipid_class"] == "PC"


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
    assert planned_jobs(snakemake(root, "--dry-run")) == {"score_spectra": 1, **AFTER_NETWORK}
    result = snakemake(root)
    assert result.returncode == 0, result.stderr[-3000:]
    assert (root / "results" / "network" / "nodes.parquet").exists()



def test_node_radius_doubles_per_decade():
    from atlas_ms.network.graph import ABSENT_RADIUS, MAX_RADIUS, MEDIAN_RADIUS, node_radius, size_reference

    reference = size_reference([1e5, 1e6, 1e7, np.nan, 0.0])
    assert reference == 1e6  # median of the detected intensities
    radii = node_radius([1e6, 1e7, 1e5, 1e12, 0.0, np.nan], reference)
    assert radii == pytest.approx([MEDIAN_RADIUS, 2 * MEDIAN_RADIUS, MEDIAN_RADIUS / 2, MAX_RADIUS,
                                   ABSENT_RADIUS, ABSENT_RADIUS])


def clustered_family(start: int = 0) -> nx.Graph:
    """
    A family shaped like a real one: dense clusters (spectra of one lipid
    class) joined by short bridges, plus a pendant chain.
    """
    rng = np.random.default_rng(1)
    graph, clusters = nx.Graph(), []
    for size in (20, 14, 8):
        nodes = list(range(start, start + size))
        graph.add_edges_from((i, j) for i in nodes for j in nodes if i < j and rng.random() < 0.5)
        clusters.append(nodes)
        start += size
    nx.add_path(graph, [clusters[0][0], start, clusters[1][0]])
    nx.add_path(graph, [clusters[1][1], start + 1, start + 2, clusters[2][0]])
    nx.add_path(graph, [clusters[0][1], start + 3, start + 4, start + 5])
    return graph


def test_layout_is_compact_and_nodes_never_overlap():
    from atlas_ms.network.graph import NODE_GAP, layout

    # A clustered family, two cliques, and singletons; nodes of random sizes.
    graph = clustered_family()
    start = max(graph) + 1
    for size in (12, 5):
        graph.add_edges_from((start + i, start + j) for i in range(size) for j in range(i + 1, size))
        start += size
    graph.add_nodes_from(range(start, start + 8))
    radius = dict(zip(graph, np.random.default_rng(0).uniform(0.125, 2.0, len(graph))))
    positions = layout(graph, radius)
    assert set(positions) == set(graph)

    # No two nodes overlap, anywhere: the gap between their edges is at least NODE_GAP.
    nodes = list(positions)
    coords = np.array([positions[n] for n in nodes])
    radii = np.array([radius[n] for n in nodes])
    dist = np.sqrt(((coords[:, None] - coords[None]) ** 2).sum(-1))
    np.fill_diagonal(dist, np.inf)
    assert (dist - radii[:, None] - radii[None, :]).min() >= NODE_GAP * 0.95

    # The bounding boxes of different families (node sizes included) do not overlap.
    families = [list(f) for f in nx.connected_components(graph) if len(f) > 1]
    boxes = [(np.array([positions[n] for n in f]), np.array([radius[n] for n in f])) for f in families]
    boxes = [((c - r[:, None]).min(0), (c + r[:, None]).max(0)) for c, r in boxes]
    for i, (low_a, high_a) in enumerate(boxes):
        for low_b, high_b in boxes[i + 1:]:
            assert (high_a < low_b).any() or (high_b < low_a).any()


def test_clustered_family_layout_is_compact():
    from atlas_ms.network.graph import MEDIAN_RADIUS, layout

    # With every node at the median size (1 unit wide), the family of dense
    # clusters takes about 4 square units per node; the first layout (scaled
    # on the closest neighbours) took about 30, most of it empty.
    graph = clustered_family()
    coords = np.array(list(layout(graph).values()))
    width, height = np.ptp(coords, axis=0) + 2 * MEDIAN_RADIUS
    assert width * height / len(graph) < 8
    assert width >= height  # longest axis horizontal
