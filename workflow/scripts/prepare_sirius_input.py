"""Rule `prepare_sirius_input`: MS2, isotope patterns and adducts for SIRIUS."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.annotation.sirius_input import run_sirius_input
from atlas_ms.logs import log_to_file

log_to_file(snakemake.log[0])
run_sirius_input(
    features_file=snakemake.input.features,
    quant_file=snakemake.input.quant,
    mgf_file=snakemake.input.mgf,
    samples=snakemake.params.samples,
    out=snakemake.output[0],
    polarity=snakemake.params.polarity,
    ppm=snakemake.params.ppm,
)
