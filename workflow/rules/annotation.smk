# Annotation. Each source writes results/annotations/<source>.parquet in the
# shared candidate format (atlas_ms.annotation.schema); `harmonize` combines
# them into confidence levels and one best annotation per feature.

from atlas_ms.annotation.lipids import rules_path
from atlas_ms.annotation.ms2query import DOWNLOAD_DONE as MS2QUERY_DOWNLOAD_DONE
from atlas_ms.annotation.ms2query import models_input as ms2query_models_input

LIBRARY_NAMES = [library["name"] for library in CFG.library_search.libraries]

# The sources the harmonization combines: always the libraries and the lipid
# rules; SIRIUS and MS2Query when switched on.
ANNOTATIONS = [RESULTS["library_annotations"], RESULTS["lipid_annotations"]]
if CFG.sirius.enabled:
    ANNOTATIONS.append("results/annotations/sirius.parquet")
if CFG.ms2query.enabled:
    ANNOTATIONS.append("results/annotations/ms2query.parquet")


# Library settings each step uses (see CLAUDE.md, rule parameters).
LIBRARY_CLEANING = ("min_library_peaks", "repair_annotations")
LIBRARY_SEARCH = ("precursor_tolerance_ppm", "fragment_tolerance_da", "min_score", "min_matched_peaks", "top_n")


rule prepare_library:
    """Read and harmonize one spectral library (once per project and library file)."""
    input:
        lambda w: CFG.library_search.library(w.library)["path"],
    output:
        library="work/annotation/libraries/{library}.parquet",
        report="work/annotation/libraries/{library}.cleaning.tsv",
    wildcard_constraints:
        library=names_regex(LIBRARY_NAMES),
    params:
        polarity=CFG.adducts.polarity,
        precursor_window_da=CFG.spectrum_qc.precursor_window_da,
        rt_unit=lambda w: CFG.library_search.library(w.library)["rt_unit"],
        # Predicted spectra skip the structure checks (see libraries.prepare_library).
        in_silico=lambda w: CFG.library_search.library(w.library)["kind"] == "in_silico",
        cleaning={key: CFG.library_search.to_dict()[key] for key in LIBRARY_CLEANING},
    log:
        "logs/prepare_library/{library}.log",
    script:
        "../scripts/prepare_library.py"


rule combine_libraries:
    """All libraries in one table, duplicates removed; cleaning report; the harmonized library as MGF."""
    input:
        libraries=expand("work/annotation/libraries/{library}.parquet", library=LIBRARY_NAMES),
        reports=expand("work/annotation/libraries/{library}.cleaning.tsv", library=LIBRARY_NAMES),
    output:
        library="work/annotation/libraries.parquet",
        summary="results/annotations/library_summary.tsv",
        cleaning="results/annotations/library_cleaning.tsv",
        # The harmonized library, to browse or reuse elsewhere.
        mgf="results/annotations/harmonized_library.mgf",
    params:
        libraries=[CFG.library_search.library(name) for name in LIBRARY_NAMES],
        remove_duplicates=CFG.library_search.remove_duplicates,
    log:
        "logs/combine_libraries.log",
    script:
        "../scripts/combine_libraries.py"


rule search_libraries:
    """Identity search of every feature MS2 spectrum in the spectral libraries."""
    input:
        mgf=RESULTS["mgf"],
        features=RESULTS["features"],
        library="work/annotation/libraries.parquet",
    output:
        RESULTS["library_annotations"],
    params:
        libraries=[CFG.library_search.library(name) for name in LIBRARY_NAMES],
        search={key: CFG.library_search.to_dict()[key] for key in LIBRARY_SEARCH},
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
        annotations=ANNOTATIONS,
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


# ---- SIRIUS and MS2Query (own conda environments, off by default) ----------

rule prepare_sirius_input:
    """MS2 spectra, MS1 isotope patterns and adducts of the features, for SIRIUS."""
    input:
        features=RESULTS["features"],
        quant=RESULTS["quant"],
        mgf=RESULTS["mgf"],
        mzml=expand("work/mzml/{sample}.mzML", sample=NAMES),
        trafo=expand("work/alignment/{sample}.trafoXML", sample=NAMES),
    output:
        "work/sirius/input.json",
    params:
        samples=NAMES,
        polarity=CFG.adducts.polarity,
        ppm=CFG.instrument.mass_error_ppm,
    log:
        "logs/prepare_sirius_input.log",
    script:
        "../scripts/prepare_sirius_input.py"


rule run_sirius:
    """SIRIUS 6: formula + ZODIAC, El Gordo, CSI:FingerID, CANOPUS (REST API)."""
    input:
        "work/sirius/input.json",
    output:
        candidates="results/annotations/sirius.parquet",
        # Kept to be opened in the SIRIUS GUI.
        project="work/sirius/project.sirius",
    # SIRIUS computes in parallel itself: nothing else runs meanwhile.
    threads: workflow.cores
    params:
        sirius=CFG.sirius.to_dict(),
        adducts=[adduct["name"] for adduct in CFG.adducts.adducts],
    log:
        "logs/run_sirius.log",
    conda:
        "../envs/sirius.yaml"
    script:
        "../scripts/run_sirius.py"


rule download_ms2query_models:
    """MS2Query's positive-mode library and models, downloaded once per machine."""
    output:
        str(MS2QUERY_DOWNLOAD_DONE),
    log:
        "logs/download_ms2query_models.log",
    conda:
        "../envs/ms2query.yaml"
    script:
        "../scripts/download_ms2query_models.py"


rule run_ms2query:
    """MS2Query analog search of every feature MS2 spectrum."""
    input:
        mgf=RESULTS["mgf"],
        models=ms2query_models_input(CFG.ms2query.models_dir),
    output:
        candidates="results/annotations/ms2query.parquet",
        csv="work/ms2query/ms2_spectra.csv",
    threads: workflow.cores
    params:
        ms2query={key: value for key, value in CFG.ms2query.to_dict().items() if key not in ("enabled", "models_dir")},
    log:
        "logs/run_ms2query.log",
    conda:
        "../envs/ms2query.yaml"
    script:
        "../scripts/run_ms2query.py"
