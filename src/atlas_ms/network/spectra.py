"""
MS2 spectra for networking: loading and quality control.

The network uses one MS2 spectrum per feature: the one exported to GNPS
(``results/gnps/ms2_spectra.mgf``, where ``SCANS`` is the feature id). Each
spectrum is cleaned, then checked with filters that do not depend on any
similarity score (see CLAUDE.md, "never filter on the score"):

* cleaning: fragments near the precursor are removed (the unfragmented
  precursor and its isotopes carry no structural information and would make
  every pair of spectra with close precursors look similar), and
  intensities are normalised to the base peak;
* too few fragments: a spectrum with only one or two peaks makes
  unreliable edges (checked here, before scoring);
* blank features: a feature that is about as intense in the blanks as in
  the samples is contamination, and its edges would only link noise
  (checked when the network is built, ``graph.py``, so that editing the
  sample types does not re-run the slow scoring).

Spectra that fail stay in the network as nodes without spectral edges, so no
feature disappears. They can still be linked by ion-identity edges.
"""

import logging
import pickle
from pathlib import Path

import pandas as pd
from matchms import Spectrum
from matchms import filtering as msfilters
from matchms.importing import load_from_mgf

from atlas_ms.config import SpectrumQCSettings

log = logging.getLogger(__name__)


def load_feature_spectra(mgf_file: str | Path) -> list[Spectrum]:
    """
    Read the exported MGF. Each spectrum gets the integer ``feature_id``
    (the MGF ``SCANS`` field); OpenMS's own FEATURE_ID is not used.
    """
    spectra = []
    for spectrum in load_from_mgf(str(mgf_file)):
        spectrum.set("feature_id", int(spectrum.get("scans")))
        spectra.append(spectrum)
    return spectra


def clean_spectrum(spectrum: Spectrum, polarity: str, settings: SpectrumQCSettings) -> Spectrum:
    """
    Standard matchms clean-up, removal of the precursor region and intensity
    normalisation. The ion mode is recorded because MS2DeepScore uses it.
    """
    spectrum = msfilters.default_filters(spectrum)
    spectrum.set("ionmode", polarity)
    spectrum = msfilters.remove_peaks_around_precursor_mz(spectrum, mz_tolerance=settings.precursor_window_da)
    return msfilters.normalize_intensities(spectrum)


def spectrum_qc(
    spectra: list[Spectrum], polarity: str, settings: SpectrumQCSettings
) -> tuple[list[Spectrum], pd.DataFrame]:
    """
    Clean all spectra and keep those with enough fragments for scoring.

    Returns the cleaned spectra that passed, and a QC table with one row per
    spectrum: feature_id, n_peaks (after cleaning), enough_peaks.
    """
    passed, rows = [], []
    for spectrum in spectra:
        cleaned = clean_spectrum(spectrum, polarity, settings)
        n_peaks = 0 if cleaned is None else len(cleaned.peaks.mz)
        enough = n_peaks >= settings.min_peaks
        rows.append({"feature_id": spectrum.get("feature_id"), "n_peaks": n_peaks, "enough_peaks": enough})
        if enough:
            passed.append(cleaned)
    qc = pd.DataFrame(rows, columns=["feature_id", "n_peaks", "enough_peaks"])
    log.info("%d spectra, %d with at least %d fragments after cleaning", len(qc), len(passed), settings.min_peaks)
    return passed, qc


def run_spectrum_qc(
    mgf_file: str | Path,
    spectra_out: str | Path,
    qc_out: str | Path,
    polarity: str,
    settings: SpectrumQCSettings,
) -> None:
    """File-level entry point: cleaned spectra (pickle, for scoring) and QC table."""
    passed, qc = spectrum_qc(load_feature_spectra(mgf_file), polarity, settings)
    with open(spectra_out, "wb") as handle:
        pickle.dump(passed, handle)
    qc.to_parquet(qc_out, index=False)


def load_spectra(path: str | Path) -> list[Spectrum]:
    """The cleaned spectra written by ``run_spectrum_qc``."""
    with open(path, "rb") as handle:
        return pickle.load(handle)
