"""Rule `annotate_run`: adduct grouping and MS2 mapping of one run."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.config import AdductSettings
from atlas_ms.logs import log_to_file
from atlas_ms.preprocessing.annotate import annotate_run

log_to_file(snakemake.log[0])
annotate_run(
    mzml=snakemake.input.mzml,
    precursors=snakemake.input.precursors,
    features_in=snakemake.input.features,
    trafo=snakemake.input.trafo,
    features_out=snakemake.output[0],
    run_name=snakemake.wildcards.sample,
    adducts=AdductSettings(**snakemake.params.adducts),
)
