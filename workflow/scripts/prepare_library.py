"""Rule `prepare_library`: read, clean and sort one spectral library."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.annotation.libraries import run_prepare_library
from atlas_ms.logs import log_to_file

log_to_file(snakemake.log[0])
run_prepare_library(
    path=snakemake.input[0],
    out=snakemake.output[0],
    polarity=snakemake.params.polarity,
    precursor_window_da=snakemake.params.precursor_window_da,
    rt_unit=snakemake.params.rt_unit,
)
