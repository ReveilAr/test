"""Rule `plan_gap_filling`: complete features and re-extraction targets."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.logs import log_to_file
from atlas_ms.preprocessing.gap_filling import plan_gap_filling

log_to_file(snakemake.log[0])
plan_gap_filling(
    consensus_file=snakemake.input[0],
    targets_out=snakemake.output.targets,
    complete_out=snakemake.output.complete,
)
