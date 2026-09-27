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
        # The feature maps that were linked: their MS2 spectra, and which
        # values gap filling re-extracted.
        features=(
            expand("work/gap_filling/{sample}.featureXML", sample=NAMES)
            if CFG.gap_filling.enabled
            else expand("work/annotated/{sample}.featureXML", sample=NAMES)
        ),
        mzml=expand("work/mzml/{sample}.mzML", sample=NAMES),
    output:
        features=RESULTS["features"],
        quant=RESULTS["quant"],
        quant_gap_filled=RESULTS["quant_gap_filled"],
        mgf=RESULTS["mgf"],
        gnps_quant=RESULTS["gnps_quant"],
        gnps_pairs=RESULTS["gnps_pairs"],
    params:
        names=NAMES,
        export=CFG.export.to_dict(),
    log:
        "logs/export.log",
    script:
        "../scripts/export.py"


rule gnps_metadata:
    """
    Sample metadata for GNPS. Its own rule, so that editing the metadata
    (sample types, ATTRIBUTE_ columns) re-runs nothing else.
    """
    input:
        "samples.tsv",
    output:
        RESULTS["gnps_metadata"],
    params:
        mzml=expand("work/mzml/{sample}.mzML", sample=NAMES),
    log:
        "logs/gnps_metadata.log",
    script:
        "../scripts/gnps_metadata.py"
