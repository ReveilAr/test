"""
Molecular network construction (the GNPS / FBMN algorithm).

From the candidate edges (``scoring.py``):

1. keep edges with score >= ``min_score`` (and, for the modified cosine, at
   least ``min_matched_peaks`` matched fragments);
2. top-K: keep an edge only if each of its two nodes is among the other's
   ``top_k`` best neighbours. This stops a node with many similar spectra
   (e.g. every PC sharing the m/z 184 fragment) from becoming a hub;
3. family size limit: while a family (connected component) has more than
   ``max_family_size`` nodes, its weakest edges are removed;
4. ion-identity edges (adducts of the same molecule, found during
   preprocessing) are added, as in Ion Identity Molecular Networking.

Nodes are the features with MS2 (the GNPS export). Each node gets:

* ``family``: its connected component, numbered by size (1 = largest);
  -1 for nodes without any edge (singletons), as in GNPS;
* ``community``: Louvain communities inside the families, finer groups of
  closely related spectra; -1 for singletons;
* ``x``, ``y``: layout coordinates, computed once here so the app shows the
  network instantly.
"""

import logging
import math

from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

from atlas_ms.config import NetworkSettings
from atlas_ms.project import read_samples

log = logging.getLogger(__name__)

EDGE_COLUMNS = ["source", "target", "edge_type", "score", "matched_peaks", "mz_delta"]


def blank_ratios(quant: pd.DataFrame, samples: pd.DataFrame) -> pd.Series:
    """
    Mean intensity in blanks / mean intensity in samples, per feature id.
    Missing values count as 0 (not detected). NaN when there are no blanks.
    """
    blanks = samples.index[samples["sample_type"] == "blank"]
    real = samples.index[samples["sample_type"] == "sample"]
    values = quant.set_index("feature_id").fillna(0.0)
    if len(blanks) == 0 or len(real) == 0:
        return pd.Series(np.nan, index=values.index)
    mean_samples = values[real].mean(axis=1)
    # A feature absent from all samples but present in blanks: ratio infinite.
    return values[blanks].mean(axis=1) / mean_samples.replace(0.0, np.nan).fillna(1e-12)


def filter_edges(candidates: pd.DataFrame, settings: NetworkSettings) -> pd.DataFrame:
    """Steps 1 and 2: score cutoffs, then top-K on both ends of each edge."""
    edges = candidates[candidates["score"] >= settings.min_score]
    if edges["matched_peaks"].notna().any():  # modified cosine
        edges = edges[edges["matched_peaks"] >= settings.min_matched_peaks]
    if edges.empty:
        return edges

    # Rank every edge among the edges of each of its two nodes (1 = best).
    ends = pd.concat([
        pd.DataFrame({"node": edges["source"], "score": edges["score"], "edge": edges.index}),
        pd.DataFrame({"node": edges["target"], "score": edges["score"], "edge": edges.index}),
    ])
    ends["rank"] = ends.groupby("node")["score"].rank(ascending=False, method="first")
    # An edge survives if it is in the top K of *both* nodes (its worst rank <= K).
    worst_rank = ends.groupby("edge")["rank"].max()
    return edges[worst_rank.loc[edges.index] <= settings.top_k]


def limit_family_size(edges: pd.DataFrame, max_size: int) -> pd.DataFrame:
    """
    Step 3: remove the weakest edges of every family larger than
    ``max_size`` until it splits into small enough families. 0 disables it.
    """
    if max_size <= 0 or edges.empty:
        return edges
    graph = nx.Graph()
    for row in edges.itertuples():
        graph.add_edge(row.source, row.target, score=row.score)
    n_removed = 0
    while True:
        too_big = [family for family in nx.connected_components(graph) if len(family) > max_size]
        if not too_big:
            break
        for family in too_big:
            family_edges = list(graph.subgraph(family).edges(data="score"))
            weakest = min(score for _, _, score in family_edges)
            doomed = [(u, v) for u, v, score in family_edges if score <= weakest]
            graph.remove_edges_from(doomed)
            n_removed += len(doomed)
    if n_removed:
        log.info("%d weak edges removed to keep families at most %d nodes", n_removed, max_size)
    kept = {tuple(sorted(edge)) for edge in graph.edges()}
    return edges[[(s, t) in kept for s, t in zip(edges["source"], edges["target"])]]


def ion_identity_edges(features: pd.DataFrame) -> pd.DataFrame:
    """
    Step 4: one edge between every two features that are adducts of the same
    molecule (``ion_partners`` from preprocessing), between features with MS2.
    """
    with_ms2 = set(features.loc[features["has_ms2"], "feature_id"])
    pairs = set()
    for row in features.dropna(subset=["ion_partners"]).itertuples():
        for partner in str(row.ion_partners).split(";"):
            partner = int(partner)
            if row.feature_id in with_ms2 and partner in with_ms2:
                pairs.add((min(row.feature_id, partner), max(row.feature_id, partner)))
    return pd.DataFrame(sorted(pairs), columns=["source", "target"]).assign(
        edge_type="ion_identity", score=np.nan, matched_peaks=np.nan
    )


def number_groups(groups: list[set]) -> dict:
    """Node -> group number, 1 for the largest group; singletons get -1."""
    numbering = {}
    ordered = sorted((g for g in groups if len(g) > 1), key=lambda g: (-len(g), min(g)))
    for number, group in enumerate(ordered, start=1):
        for node in group:
            numbering[node] = number
    return numbering


def layout(graph: nx.Graph) -> dict:
    """
    Node -> (x, y). Each family gets its own force-directed layout in a box
    whose side grows with the square root of its size. The boxes are placed
    row by row, largest first, with the singletons at the end, like
    Cytoscape's grid of components.
    """
    families = sorted(nx.connected_components(graph), key=lambda g: (-len(g), min(g)))
    total_area = sum(len(family) for family in families)
    row_width = max(math.sqrt(total_area) * 2.0, math.sqrt(len(families[0])) * 2.0) if families else 1.0
    positions, x, y, row_height = {}, 0.0, 0.0, 0.0
    for family in families:
        side = math.sqrt(len(family)) * 2.0 if len(family) > 1 else 1.0
        if x + side > row_width and x > 0:  # start a new row
            x, y, row_height = 0.0, y - row_height - 1.0, 0.0
        if len(family) == 1:
            local = {next(iter(family)): np.array([0.5, 0.5])}
        else:
            spring = nx.spring_layout(graph.subgraph(family), weight="score", seed=0)
            coords = np.array(list(spring.values()))
            span = np.ptp(coords, axis=0)
            span[span == 0] = 1.0
            local = dict(zip(spring, (coords - coords.min(axis=0)) / span))
        for node, (u, v) in local.items():
            positions[node] = (x + u * side, y - v * side)
        x += side + 1.0
        row_height = max(row_height, side)
    return positions


def build_network(
    candidates: pd.DataFrame,
    features: pd.DataFrame,
    settings: NetworkSettings,
    excluded: set = frozenset(),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    The network: (nodes, edges) tables.

    ``features`` is ``results/features.parquet``; its features with MS2 are
    the nodes. Features in ``excluded`` (blank features) get no spectral
    edges; they are removed before the top-K step, so they take no place
    among anyone's best neighbours. Edge ``mz_delta`` is the precursor m/z
    difference (target - source), which hints at the chemical modification
    between two nodes.
    """
    candidates = candidates[~candidates["source"].isin(excluded) & ~candidates["target"].isin(excluded)]
    spectral = limit_family_size(filter_edges(candidates, settings), settings.max_family_size)
    spectral = spectral.assign(edge_type="spectral")
    edges = pd.concat([spectral, ion_identity_edges(features)], ignore_index=True)
    mz = features.set_index("feature_id")["mz"]
    edges["mz_delta"] = mz.loc[edges["target"]].to_numpy() - mz.loc[edges["source"]].to_numpy()
    edges = edges[EDGE_COLUMNS]

    node_ids = features.loc[features["has_ms2"], "feature_id"].tolist()
    graph = nx.Graph()
    graph.add_nodes_from(node_ids)
    for row in edges.itertuples():
        # Ion-identity edges count as strong links for communities and layout.
        weight = 1.0 if row.edge_type == "ion_identity" else row.score
        graph.add_edge(row.source, row.target, score=weight)

    family = number_groups(list(nx.connected_components(graph)))
    communities = nx.community.louvain_communities(
        graph, weight="score", resolution=settings.louvain_resolution, seed=0
    )
    community = number_groups(communities)
    positions = layout(graph)
    nodes = pd.DataFrame({
        "feature_id": node_ids,
        "family": [family.get(n, -1) for n in node_ids],
        "community": [community.get(n, -1) for n in node_ids],
        "degree": [graph.degree(n) for n in node_ids],
        "x": [positions[n][0] for n in node_ids],
        "y": [positions[n][1] for n in node_ids],
    })
    n_families = len({f for f in family.values()})
    log.info(
        "Network: %d nodes, %d spectral + %d ion-identity edges, %d families, %d singletons",
        len(nodes), int((edges["edge_type"] == "spectral").sum()),
        int((edges["edge_type"] == "ion_identity").sum()), n_families, int((nodes["family"] == -1).sum()),
    )
    return nodes, edges


def run_network(
    candidates_file: str | Path,
    features_file: str | Path,
    qc_file: str | Path,
    quant_file: str | Path,
    samples_file: str | Path,
    nodes_out: str | Path,
    edges_out: str | Path,
    graphml_out: str | Path,
    max_blank_ratio: float,
    settings: NetworkSettings,
) -> None:
    """
    File-level entry point: blank check, network construction, and the
    nodes / edges tables and GraphML file. The nodes table also carries the
    spectrum QC: n_peaks, blank_ratio, and ``qc`` (why a node has no spectral
    edges; empty when it may have them).
    """
    features = pd.read_parquet(features_file)
    ratios = blank_ratios(pd.read_parquet(quant_file), read_samples(samples_file))
    blank = set(ratios.index[ratios > max_blank_ratio])  # empty without blanks (NaN ratios)
    nodes, edges = build_network(pd.read_parquet(candidates_file), features, settings, excluded=blank)

    qc = pd.read_parquet(qc_file)
    nodes = nodes.merge(qc, on="feature_id", how="left")
    nodes["blank_ratio"] = ratios.reindex(nodes["feature_id"]).to_numpy()
    nodes["qc"] = np.where(~nodes["enough_peaks"].fillna(False).astype(bool), "too few fragments",
                           np.where(nodes["feature_id"].isin(blank), "present in blanks", ""))
    nodes = nodes.drop(columns="enough_peaks")
    log.info("%d nodes without spectral edges because of QC (%d blank features)",
             int((nodes["qc"] != "").sum()), len(blank & set(nodes["feature_id"])))

    nodes.to_parquet(nodes_out, index=False)
    edges.to_parquet(edges_out, index=False)
    to_graphml(nodes, edges, features, graphml_out)


def to_graphml(nodes: pd.DataFrame, edges: pd.DataFrame, features: pd.DataFrame, path) -> None:
    """
    Write the network as GraphML (for Cytoscape), with the feature table's
    columns as node attributes. GraphML has no missing values: they become
    empty strings (text) or NaN (numbers).
    """
    table = nodes.merge(features, on="feature_id", how="left")
    graph = nx.Graph()
    for record in table.to_dict("records"):
        attributes = {key: _graphml_value(value) for key, value in record.items()}
        graph.add_node(int(record["feature_id"]), **attributes)
    for record in edges.to_dict("records"):
        attributes = {key: _graphml_value(value) for key, value in record.items() if key not in ("source", "target")}
        graph.add_edge(int(record["source"]), int(record["target"]), **attributes)
    nx.write_graphml(graph, path)


def _graphml_value(value):
    """Convert a table value into a type GraphML accepts."""
    if value is None:
        return ""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value)
    return str(value)
