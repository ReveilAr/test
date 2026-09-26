"""Rule `run_sirius` (sirius conda environment): annotation with SIRIUS 6."""
# `snakemake` is provided by Snakemake's script directive.
import sys
from pathlib import Path

# atlas_ms is not installed in this environment: use the clone's source
# (Snakemake sets __file__ to this script's real path). Only light modules
# are imported from it here: no pyOpenMS, no matchms.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from atlas_ms.annotation.sirius import run_sirius_job  # noqa: E402
from atlas_ms.logs import log_to_file  # noqa: E402

log_to_file(snakemake.log[0])
run_sirius_job(
    input_file=snakemake.input[0],
    out=snakemake.output.candidates,
    project_file=snakemake.output.project,
    settings=snakemake.params.sirius,
    fallback_adducts=snakemake.params.adducts,
)
