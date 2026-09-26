"""Rule `download_ms2query_models` (ms2query conda environment)."""
# `snakemake` is provided by Snakemake's script directive.
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # atlas_ms is not installed here

from atlas_ms.annotation.ms2query import download_models  # noqa: E402
from atlas_ms.logs import log_to_file  # noqa: E402

log_to_file(snakemake.log[0])
download_models(snakemake.output[0])
