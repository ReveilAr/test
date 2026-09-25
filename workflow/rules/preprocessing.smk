# Feature finding, alignment, adduct grouping, MS2 mapping and linking
# (UmetaFlow preprocessing, ported to pyOpenMS).


rule find_features:
    """Untargeted feature detection in one run + MS2 precursor correction."""
    input:
        mzml="work/mzml/{sample}.mzML",
    output:
        features="work/features/{sample}.featureXML",
        precursors="work/features/{sample}.precursors.tsv",
    params:
        instrument=CFG.instrument.to_dict(),
        feature_finding=CFG.feature_finding.to_dict(),
    log:
        "logs/find_features/{sample}.log",
    script:
        "../scripts/find_features.py"


rule align:
    """Retention-time alignment: one transformation per run."""
    input:
        expand("work/features/{sample}.featureXML", sample=NAMES),
    output:
        expand("work/alignment/{sample}.trafoXML", sample=NAMES),
    params:
        alignment=CFG.alignment.to_dict(),
    log:
        "logs/align.log",
    script:
        "../scripts/align.py"


rule annotate_run:
    """Adduct grouping and MS2 mapping of one run, on the aligned RT axis."""
    input:
        mzml="work/mzml/{sample}.mzML",
        precursors="work/features/{sample}.precursors.tsv",
        features="work/features/{sample}.featureXML",
        trafo="work/alignment/{sample}.trafoXML",
    output:
        "work/annotated/{sample}.featureXML",
    params:
        adducts=CFG.adducts.to_dict(),
    log:
        "logs/annotate_run/{sample}.log",
    script:
        "../scripts/annotate_run.py"


rule link:
    """Group corresponding features of all runs into consensus features."""
    input:
        features=expand("work/annotated/{sample}.featureXML", sample=NAMES),
        mzml=expand("work/mzml/{sample}.mzML", sample=NAMES),
    output:
        "work/consensus/linked.consensusXML",
    params:
        names=NAMES,
        linking=CFG.linking.to_dict(),
    log:
        "logs/link.log",
    script:
        "../scripts/link.py"
