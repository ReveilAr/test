"""Rule `build_network`: molecular network from the candidate edges."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.config import NetworkSettings
from atlas_ms.logs import log_to_file
from atlas_ms.network.graph import run_network

log_to_file(snakemake.log[0])
run_network(
    candidates_file=snakemake.input.candidates,
    features_file=snakemake.input.features,
    qc_file=snakemake.input.qc,
    quant_file=snakemake.input.quant,
    samples_file=snakemake.input.samples,
    nodes_out=snakemake.output.nodes,
    edges_out=snakemake.output.edges,
    graphml_out=snakemake.output.graphml,
    max_blank_ratio=snakemake.params.max_blank_ratio,
    settings=NetworkSettings(**snakemake.params.network),
)
