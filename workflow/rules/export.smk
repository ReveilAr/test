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
        # Metadata edits (sample types, ATTRIBUTE_ columns) only re-run this rule.
        samples="samples.tsv",
    output:
        features=RESULTS["features"],
        quant=RESULTS["quant"],
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
