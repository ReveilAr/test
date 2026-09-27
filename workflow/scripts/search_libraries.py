"""Rule `search_libraries`: identity search of the feature spectra in the libraries."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.annotation.libraries import run_search_libraries
from atlas_ms.config import LibrarySearchSettings
from atlas_ms.logs import log_to_file

log_to_file(snakemake.log[0])
run_search_libraries(
    mgf_file=snakemake.input.mgf,
    features_file=snakemake.input.features,
    library_file=snakemake.input.library,
    entries=snakemake.params.libraries,
    out=snakemake.output[0],
    precursor_window_da=snakemake.params.precursor_window_da,
    settings=LibrarySearchSettings(**snakemake.params.search),
)
