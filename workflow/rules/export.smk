# Final feature tables (Parquet) and GNPS/FBMN export.


rule export:
    """Feature and quantification tables + GNPS FBMN files."""
    input:
        # With gap filling off, the first linking result is exported directly.
        consensus=(
            "work/consensus/gap_filled.consensusXML"
            if CFG.gap_filling.enabled
            else "work/consensus/linked.consensusXML"
        ),
        mzml=expand("work/mzml/{sample}.mzML", sample=NAMES),
        # Gap-filled feature maps: tell which values were re-extracted.
        gap_filled=(
            expand("work/gap_filling/{sample}.featureXML", sample=NAMES)
            if CFG.gap_filling.enabled
            else []
        ),
        # Metadata edits (sample types, ATTRIBUTE_ columns) only re-run this rule.
        samples="samples.tsv",
    output:
        features=RESULTS["features"],
        quant=RESULTS["quant"],
        quant_gap_filled=RESULTS["quant_gap_filled"],
        mgf=RESULTS["mgf"],
        gnps_quant=RESULTS["gnps_quant"],
        gnps_metadata=RESULTS["gnps_metadata"],
        gnps_pairs=RESULTS["gnps_pairs"],
        gnps_consensus="work/export/gnps.consensusXML",
    params:
        export=CFG.export.to_dict(),
    log:
        "logs/export.log",
    script:
        "../scripts/export.py"
