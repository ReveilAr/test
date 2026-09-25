"""Rule `align`: retention-time alignment of all runs."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.config import AlignmentSettings
from atlas_ms.logs import log_to_file
from atlas_ms.preprocessing.alignment import align_runs

log_to_file(snakemake.log[0])
align_runs(
    feature_files=list(snakemake.input),
    trafo_files=list(snakemake.output),
    settings=AlignmentSettings(**snakemake.params.alignment),
)
