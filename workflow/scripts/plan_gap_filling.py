"""Rule `plan_gap_filling`: re-extraction targets and feature membership."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.logs import log_to_file
from atlas_ms.preprocessing.gap_filling import plan_gap_filling

log_to_file(snakemake.log[0])
plan_gap_filling(
    consensus_file=snakemake.input[0],
    targets_out=snakemake.output.targets,
    members_out=snakemake.output.members,
)
