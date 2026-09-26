# Annotation. Each source writes results/annotations/<source>.parquet in the
# shared candidate format (atlas_ms.annotation.schema); `harmonize` combines
# them into confidence levels and one best annotation per feature.

from atlas_ms.annotation.lipids import rules_path

LIBRARY_NAMES = [library["name"] for library in CFG.library_search.libraries]


rule prepare_library:
    """Read, clean and sort one spectral library (once per project and library file)."""
    input:
        lambda w: CFG.library_search.library(w.library)["path"],
    output:
        "work/annotation/libraries/{library}.parquet",
    wildcard_constraints:
        library=names_regex(LIBRARY_NAMES),
    params:
        polarity=CFG.adducts.polarity,
        precursor_window_da=CFG.spectrum_qc.precursor_window_da,
        rt_unit=lambda w: CFG.library_search.library(w.library)["rt_unit"],
    log:
        "logs/prepare_library/{library}.log",
    script:
        "../scripts/prepare_library.py"


rule search_libraries:
    """Identity search of every feature MS2 spectrum in the spectral libraries."""
    input:
        mgf=RESULTS["mgf"],
        features=RESULTS["features"],
        libraries=expand("work/annotation/libraries/{library}.parquet", library=LIBRARY_NAMES),
    output:
        RESULTS["library_annotations"],
    params:
        libraries=[CFG.library_search.library(name) for name in LIBRARY_NAMES],
        search={key: value for key, value in CFG.library_search.to_dict().items() if key != "libraries"},
        precursor_window_da=CFG.spectrum_qc.precursor_window_da,
    log:
        "logs/search_libraries.log",
    script:
        "../scripts/search_libraries.py"


rule annotate_lipids:
    """Rule-based lipid annotation (class fragments and losses, species, chains)."""
    input:
        mgf=RESULTS["mgf"],
        features=RESULTS["features"],
        rules=rules_path(CFG.lipids.rules_file),
    output:
        RESULTS["lipid_annotations"],
    params:
        lipids={key: value for key, value in CFG.lipids.to_dict().items() if key != "rules_file"},
    log:
        "logs/annotate_lipids.log",
    script:
        "../scripts/annotate_lipids.py"


rule harmonize:
    """Confidence levels, flags, best annotation per feature, family class consensus."""
    input:
        annotations=[RESULTS["library_annotations"], RESULTS["lipid_annotations"]],
        features=RESULTS["features"],
        # Families and blank flags.
        nodes=RESULTS["nodes"],
    output:
        candidates=RESULTS["candidates"],
        best=RESULTS["best"],
    params:
        harmonization=CFG.harmonization.to_dict(),
    log:
        "logs/harmonize.log",
    script:
        "../scripts/harmonize.py"


rule export_graphml:
    """The network for Cytoscape, with the best annotations as node attributes."""
    input:
        nodes=RESULTS["nodes"],
        edges=RESULTS["edges"],
        features=RESULTS["features"],
        best=RESULTS["best"],
    output:
        RESULTS["graphml"],
    log:
        "logs/export_graphml.log",
    script:
        "../scripts/export_graphml.py"
