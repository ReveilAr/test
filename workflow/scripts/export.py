"""Rule `export`: feature tables and GNPS/FBMN files."""
# `snakemake` is provided by Snakemake's script directive.
from pathlib import Path

from atlas_ms.config import ExportSettings
from atlas_ms.logs import log_to_file
from atlas_ms.preprocessing.export import export_results

log_to_file(snakemake.log[0])
export_results(
    consensus_file=snakemake.input.consensus,
    feature_files=list(snakemake.input.features),
    mzml_files=list(snakemake.input.mzml),
    run_names=snakemake.params.names,
    features_out=snakemake.output.features,
    quant_out=snakemake.output.quant,
    gap_filled_out=snakemake.output.quant_gap_filled,
    gnps_dir=Path(snakemake.output.mgf).parent,
    settings=ExportSettings(**snakemake.params.export),
)
