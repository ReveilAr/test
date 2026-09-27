# Molecular network: spectrum QC -> pairwise scoring -> network construction.
# Scoring is the slow step and depends only on the spectra and the scoring
# parameters; the network cutoffs, blank check and layout are in the cheap
# last step.

from atlas_ms.network.model_files import DEFAULT_MODEL_PATH, ms2deepscore_model_path


def scoring_params():
    """Only the scoring values the chosen score uses (see CLAUDE.md, rule parameters)."""
    values = CFG.scoring.to_dict()
    if CFG.scoring.score == "modified_cosine":
        values.pop("ms2deepscore_model")
    else:
        values.pop("fragment_tolerance_da")
    return values


rule spectrum_qc:
    """Clean the feature MS2 spectra and keep those with enough fragments."""
    input:
        RESULTS["mgf"],
    output:
        spectra="work/network/spectra.pickle",
        qc="work/network/spectrum_qc.parquet",
    params:
        polarity=CFG.adducts.polarity,
        spectrum_qc={
            "precursor_window_da": CFG.spectrum_qc.precursor_window_da,
            "min_peaks": CFG.spectrum_qc.min_peaks,
        },
    log:
        "logs/spectrum_qc.log",
    script:
        "../scripts/spectrum_qc.py"


rule score_spectra:
    """All-against-all similarity -> pool of candidate edges."""
    input:
        spectra="work/network/spectra.pickle",
        model=(
            str(ms2deepscore_model_path(CFG.scoring))
            if CFG.scoring.score == "ms2deepscore"
            else []
        ),
    output:
        "work/network/candidates.parquet",
    # MS2DeepScore (PyTorch) computes in parallel: all cores. The modified
    # cosine runs on one core.
    threads: workflow.cores if CFG.scoring.score == "ms2deepscore" else 1
    params:
        scoring=scoring_params(),
    log:
        "logs/score_spectra.log",
    script:
        "../scripts/score_spectra.py"


rule build_network:
    """Blank check, GNPS-style network, families, communities, layout."""
    input:
        candidates="work/network/candidates.parquet",
        features=RESULTS["features"],
        qc="work/network/spectrum_qc.parquet",
        quant=RESULTS["quant"],
        # Sample types decide which features are blank signal.
        samples="samples.tsv",
    output:
        nodes=RESULTS["nodes"],
        edges=RESULTS["edges"],
    params:
        max_blank_ratio=CFG.spectrum_qc.max_blank_ratio,
        network=CFG.network.to_dict(),
    log:
        "logs/build_network.log",
    script:
        "../scripts/build_network.py"


rule download_ms2deepscore_model:
    """The pretrained MS2DeepScore model, downloaded once per machine."""
    output:
        str(DEFAULT_MODEL_PATH),
    log:
        "logs/download_ms2deepscore_model.log",
    script:
        "../scripts/download_ms2deepscore_model.py"
