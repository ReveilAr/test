"""Rule `export_graphml`: the annotated network as GraphML (Cytoscape)."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.logs import log_to_file
from atlas_ms.network.graph import run_export_graphml

log_to_file(snakemake.log[0])
run_export_graphml(
    nodes_file=snakemake.input.nodes,
    edges_file=snakemake.input.edges,
    features_file=snakemake.input.features,
    best_file=snakemake.input.best,
    out=snakemake.output[0],
)
