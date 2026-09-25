"""Rule `spectrum_qc`: clean the feature MS2 spectra before scoring."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.config import SpectrumQCSettings
from atlas_ms.logs import log_to_file
from atlas_ms.network.spectra import run_spectrum_qc

log_to_file(snakemake.log[0])
run_spectrum_qc(
    mgf_file=snakemake.input[0],
    spectra_out=snakemake.output.spectra,
    qc_out=snakemake.output.qc,
    polarity=snakemake.params.polarity,
    settings=SpectrumQCSettings(**snakemake.params.spectrum_qc),
)
