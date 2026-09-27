"""Rule `score_spectra`: pairwise spectral similarity -> candidate edges."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.config import ScoringSettings
from atlas_ms.logs import log_to_file
from atlas_ms.network.scoring import run_scoring

log_to_file(snakemake.log[0])
run_scoring(
    spectra_file=snakemake.input.spectra,
    candidates_out=snakemake.output[0],
    settings=ScoringSettings(**snakemake.params.scoring),
    threads=snakemake.threads,
)
