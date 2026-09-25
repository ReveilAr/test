"""
LC-MS preprocessing with pyOpenMS: raw data -> aligned feature table + GNPS export.

A port of UmetaFlow's Snakemake rules (which call the OpenMS command-line
tools) to pyOpenMS, in the same order:

    features     untargeted feature finding in each run
    alignment    retention-time alignment across runs
    annotate     adduct grouping + MS2 mapping per run
    linking      grouping of features across runs
    gap_filling  targeted re-extraction of missing values
    export       feature/quantification tables + GNPS FBMN files
"""
