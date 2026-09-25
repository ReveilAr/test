"""Rule `find_features`: untargeted feature detection in one run."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.config import FeatureFindingSettings, InstrumentSettings
from atlas_ms.logs import log_to_file
from atlas_ms.preprocessing.features import find_features

log_to_file(snakemake.log[0])
find_features(
    mzml=snakemake.input.mzml,
    features_out=snakemake.output.features,
    precursors_out=snakemake.output.precursors,
    instrument=InstrumentSettings(**snakemake.params.instrument),
    settings=FeatureFindingSettings(**snakemake.params.feature_finding),
)
