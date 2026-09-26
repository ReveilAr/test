"""Rule `annotate_lipids`: rule-based lipid annotation."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.annotation.lipids import run_annotate_lipids
from atlas_ms.config import LipidSettings
from atlas_ms.logs import log_to_file

log_to_file(snakemake.log[0])
run_annotate_lipids(
    mgf_file=snakemake.input.mgf,
    features_file=snakemake.input.features,
    rules_file=snakemake.input.rules,
    out=snakemake.output[0],
    settings=LipidSettings(**snakemake.params.lipids),
)
