"""Rule `export`: feature tables and GNPS/FBMN files."""
# `snakemake` is provided by Snakemake's script directive.
from pathlib import Path

from atlas_ms.config import ExportSettings
from atlas_ms.logs import log_to_file
from atlas_ms.preprocessing.export import export_results
from atlas_ms.project import read_samples

log_to_file(snakemake.log[0])
export_results(
    consensus_file=snakemake.input.consensus,
    mzml_files=list(snakemake.input.mzml),
    samples=read_samples(snakemake.input.samples),
    features_out=snakemake.output.features,
    quant_out=snakemake.output.quant,
    gap_filled_out=snakemake.output.quant_gap_filled,
    gnps_dir=Path(snakemake.output.mgf).parent,
    gnps_consensus_out=snakemake.output.gnps_consensus,
    settings=ExportSettings(**snakemake.params.export),
    gap_filled_files=list(snakemake.input.gap_filled),
)
