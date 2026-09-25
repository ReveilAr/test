"""Rule `download_ms2deepscore_model`: fetch the pretrained model once."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.logs import log_to_file
from atlas_ms.network.model_files import download_model

log_to_file(snakemake.log[0])
download_model(snakemake.output[0])
