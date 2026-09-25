"""Network construction steps (pure logic) and the network of the synthetic study."""

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from atlas_ms.config import NetworkSettings, ScoringSettings
from atlas_ms.network.graph import filter_edges, limit_family_size, number_groups
from atlas_ms.network.scoring import candidate_edges
from test_workflow import feature_of


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
    """MS2DeepScore path, with a tiny untrained model (the real one is downloaded on first use)."""
    from ms2deepscore.models import SiameseSpectralModel
    from ms2deepscore.SettingsMS2Deepscore import SettingsMS2Deepscore

    from atlas_ms.network.scoring import score_spectra
    from atlas_ms.network.spectra import load_spectra

    model_file = tmp_path / "tiny.pt"
    SiameseSpectralModel(SettingsMS2Deepscore(base_dims=(32,), embedding_dim=8)).save(str(model_file))
    spectra = load_spectra(processed_project.root / "work" / "network" / "spectra.pickle")
    settings = ScoringSettings(score="ms2deepscore", ms2deepscore_model=str(model_file))
    pool = score_spectra(spectra, settings)
    # An untrained model scores everything alike: every pair is a candidate.
    n = len(spectra)
    assert len(pool) == n * (n - 1) // 2
    assert pool["matched_peaks"].isna().all()
    assert pool["score"].between(-1, 1).all()


def test_config_rejects_inconsistent_network_cutoffs():
    from atlas_ms.config import ProjectConfig

    with pytest.raises(ValueError, match="min_candidate_score"):
        ProjectConfig.from_dict({"network": {"min_score": 0.2}})
    with pytest.raises(ValueError, match="candidates_per_spectrum"):
        ProjectConfig.from_dict({"network": {"top_k": 80}})
