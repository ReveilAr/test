"""Rule `harmonize`: combine the annotation sources."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.annotation.harmonize import run_harmonize
from atlas_ms.config import HarmonizationSettings
from atlas_ms.logs import log_to_file

log_to_file(snakemake.log[0])
run_harmonize(
    annotation_files=list(snakemake.input.annotations),
    features_file=snakemake.input.features,
    nodes_file=snakemake.input.nodes,
    candidates_out=snakemake.output.candidates,
    best_out=snakemake.output.best,
    settings=HarmonizationSettings(**snakemake.params.harmonization),
)
