"""Rules `link` and `link_gap_filled`: group features across runs."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.config import LinkingSettings
from atlas_ms.logs import log_to_file
from atlas_ms.preprocessing.linking import link_runs

log_to_file(snakemake.log[0])
link_runs(
    feature_files=list(snakemake.input.features),
    run_names=snakemake.params.names,
    mzml_files=list(snakemake.input.mzml),
    consensus_out=snakemake.output[0],
    settings=LinkingSettings(**snakemake.params.linking),
)
