"""Rule `gnps_metadata`: sample metadata in GNPS format."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.logs import log_to_file
from atlas_ms.preprocessing.export import gnps_metadata
from atlas_ms.project import read_samples

log_to_file(snakemake.log[0])
gnps_metadata(read_samples(snakemake.input[0]), snakemake.params.mzml).to_csv(
    snakemake.output[0], sep="\t", index=False
)
