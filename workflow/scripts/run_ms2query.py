"""Rule `run_ms2query` (ms2query conda environment): analog search."""
# `snakemake` is provided by Snakemake's script directive.
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # atlas_ms is not installed here

from atlas_ms.annotation.ms2query import run_ms2query  # noqa: E402
from atlas_ms.logs import log_to_file  # noqa: E402

log_to_file(snakemake.log[0])
run_ms2query(
    mgf_file=snakemake.input.mgf,
    models_dir=snakemake.input.models,
    csv_out=snakemake.output.csv,
    out=snakemake.output.candidates,
    settings=snakemake.params.ms2query,
    threads=snakemake.threads,
)
