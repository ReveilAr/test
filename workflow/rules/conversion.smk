# Raw data -> centroided mzML in work/mzml/<sample>.mzML.


rule convert_thermo:
    """
    Thermo .raw -> indexed mzML with ThermoRawFileParser. Peak picking
    (centroiding) by the Thermo library is on by default.
    """
    input:
        lambda wc: SAMPLES.loc[wc.sample, "file"],
    output:
        "work/mzml/{sample}.mzML",
    wildcard_constraints:
        sample=THERMO,
    log:
        "logs/convert/{sample}.log",
    conda:
        "../envs/conversion.yaml"
    shell:
        "thermorawfileparser --input={input:q} --output={output:q} --format=2 > {log:q} 2>&1"


rule link_mzml:
    """
    An mzML input is used as it is, through a symbolic link: no copy. It must
    already be centroided (checked by feature finding).
    """
    input:
        lambda wc: SAMPLES.loc[wc.sample, "file"],
    output:
        "work/mzml/{sample}.mzML",
    wildcard_constraints:
        sample=MZML,
    shell:
        "ln -sf {input:q} {output:q}"
