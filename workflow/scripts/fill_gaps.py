"""Rule `fill_gaps`: targeted re-extraction of missing features in one run."""
# `snakemake` is provided by Snakemake's script directive.
import os

# OpenMP threads for pyOpenMS: the cores Snakemake gave this step (set
# before pyOpenMS is imported, which starts OpenMP).
os.environ["OMP_NUM_THREADS"] = str(snakemake.threads)

from atlas_ms.config import AdductSettings, GapFillingSettings, InstrumentSettings  # noqa: E402
from atlas_ms.logs import log_to_file  # noqa: E402
from atlas_ms.preprocessing.gap_filling import fill_gaps  # noqa: E402

log_to_file(snakemake.log[0])
fill_gaps(
    mzml=snakemake.input.mzml,
    precursors=snakemake.input.precursors,
    features_in=snakemake.input.features,
    trafo=snakemake.input.trafo,
    targets_file=snakemake.input.targets,
    members_file=snakemake.input.members,
    map_index=snakemake.params.map_index,
    features_out=snakemake.output[0],
    run_name=snakemake.wildcards.sample,
    instrument=InstrumentSettings(**snakemake.params.instrument),
    adducts=AdductSettings(**snakemake.params.adducts),
    settings=GapFillingSettings(**snakemake.params.gap_filling),
)
