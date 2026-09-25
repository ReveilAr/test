# Gap filling: re-extract features missing in some runs, then link again.
# Only used when gap_filling.enabled is true (see rules/export.smk).


rule plan_gap_filling:
    """Re-extraction targets, and which consensus feature every feature belongs to."""
    input:
        "work/consensus/linked.consensusXML",
    output:
        targets="work/gap_filling/targets.tsv",
        members="work/gap_filling/members.tsv",
    log:
        "logs/plan_gap_filling.log",
    script:
        "../scripts/plan_gap_filling.py"


rule fill_gaps:
    """Targeted re-extraction in one run, merged with the features to keep."""
    input:
        mzml="work/mzml/{sample}.mzML",
        precursors="work/features/{sample}.precursors.tsv",
        features="work/features/{sample}.featureXML",
        trafo="work/alignment/{sample}.trafoXML",
        targets="work/gap_filling/targets.tsv",
        members="work/gap_filling/members.tsv",
    output:
        "work/gap_filling/{sample}.featureXML",
    params:
        map_index=lambda wc: NAMES.index(wc.sample),
        instrument=CFG.instrument.to_dict(),
        adducts=CFG.adducts.to_dict(),
        gap_filling=CFG.gap_filling.to_dict(),
    log:
        "logs/fill_gaps/{sample}.log",
    script:
        "../scripts/fill_gaps.py"


rule link_gap_filled:
    """Link the gap-filled runs."""
    input:
        features=expand("work/gap_filling/{sample}.featureXML", sample=NAMES),
        mzml=expand("work/mzml/{sample}.mzML", sample=NAMES),
    output:
        "work/consensus/gap_filled.consensusXML",
    params:
        names=NAMES,
        linking=CFG.linking.to_dict(),
    log:
        "logs/link_gap_filled.log",
    script:
        "../scripts/link.py"
