"""
Rule `annotate_run`: adduct grouping and MS2 mapping of one run (MS2
spectra stored as spectrum indices on the features: annotate.attach_ms2).
"""
# `snakemake` is provided by Snakemake's script directive.
import os

# OpenMP threads for pyOpenMS: the cores Snakemake gave this step (set
# before pyOpenMS is imported, which starts OpenMP).
os.environ["OMP_NUM_THREADS"] = str(snakemake.threads)

from atlas_ms.config import AdductSettings  # noqa: E402
from atlas_ms.logs import log_to_file  # noqa: E402
from atlas_ms.preprocessing.annotate import annotate_run  # noqa: E402

log_to_file(snakemake.log[0])
annotate_run(
    mzml=snakemake.input.mzml,
    precursors=snakemake.input.precursors,
    features_in=snakemake.input.features,
    trafo=snakemake.input.trafo,
    features_out=snakemake.output[0],
    run_name=snakemake.wildcards.sample,
    adducts=AdductSettings(**snakemake.params.adducts),
)
