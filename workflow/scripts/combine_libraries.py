"""Rule `combine_libraries`: one table of all libraries, duplicates removed."""
# `snakemake` is provided by Snakemake's script directive.
from atlas_ms.annotation.libraries import run_combine_libraries
from atlas_ms.logs import log_to_file

log_to_file(snakemake.log[0])
run_combine_libraries(
    library_files=list(snakemake.input.libraries),
    report_files=list(snakemake.input.reports),
    entries=snakemake.params.libraries,
    out=snakemake.output.library,
    summary_out=snakemake.output.summary,
    cleaning_out=snakemake.output.cleaning,
    remove_duplicates=snakemake.params.remove_duplicates,
)
