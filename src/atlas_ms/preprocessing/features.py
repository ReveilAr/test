"""
Untargeted feature detection in one run.

This is UmetaFlow's `precursorcorrection_peak`, `preprocess` and
`precursorcorrection_feature` rules (OpenMS command-line tools), ported to
pyOpenMS:

1. Precursor correction to MS1 peaks. The instrument reports the centre of
   the isolation window as the MS2 precursor m/z. It is moved to the most
   intense MS1 centroid nearby.
2. FeatureFinderMetabo, in three steps:
   * MassTraceDetection follows each m/z through consecutive MS1 scans
     (mass traces);
   * ElutionPeakDetection splits traces into chromatographic peaks;
   * FeatureFindingMetabo groups the isotope traces (M, M+1, M+2...) of one
     compound into a feature with m/z, RT, charge and intensity.
3. Precursor correction to features. Each MS2 precursor is moved to the
   monoisotopic m/z of the feature it falls in. This fixes MS2 spectra
   triggered on an isotope peak and lets them be matched to the right
   feature later.
"""

import logging
from pathlib import Path

import pyopenms as oms

from atlas_ms.config import FeatureFindingSettings, InstrumentSettings
from atlas_ms.preprocessing.msdata import (
    check_centroided,
    load_experiment,
    ms1_experiment,
    precursor_table,
    set_parameters,
    store_feature_map,
)

log = logging.getLogger(__name__)


def detect_features(
    experiment: oms.MSExperiment,
    instrument: InstrumentSettings,
    settings: FeatureFindingSettings,
) -> oms.FeatureMap:
    """Run the three FeatureFinderMetabo steps on the MS1 scans of ``experiment``."""
    ms1 = ms1_experiment(experiment, min_intensity=instrument.noise_threshold)

    # 1) Mass traces: an m/z followed over consecutive scans within mass_error_ppm.
    mass_traces = []
    detector = oms.MassTraceDetection()
    set_parameters(detector, {
        "mass_error_ppm": instrument.mass_error_ppm,
        "noise_threshold_int": instrument.noise_threshold,
        "min_trace_length": instrument.min_trace_length_s,
    })
    detector.run(ms1, mass_traces, 0)  # 0 = no limit on the number of traces

    # 2) Split traces at chromatographic minima into single elution peaks.
    elution_peaks = []
    peak_detector = oms.ElutionPeakDetection()
    set_parameters(peak_detector, {"chrom_fwhm": instrument.chrom_fwhm_s, "width_filtering": "fixed"})
    peak_detector.detectPeaks(mass_traces, elution_peaks)

    # 3) Group co-eluting isotope traces into features. Convex hulls are kept:
    #    MS2 spectra are later mapped to the feature whose hull contains them.
    features = oms.FeatureMap()
    finder = oms.FeatureFindingMetabo()
    set_parameters(finder, {
        "chrom_fwhm": instrument.chrom_fwhm_s,
        "isotope_filtering_model": settings.isotope_filtering_model,
        "remove_single_traces": "true" if settings.remove_single_traces else "false",
        "report_convex_hulls": "true",
    })
    finder.run(elution_peaks, features, [])
    features.setUniqueIds()  # stable ids, used to follow features across steps

    log.info("%d mass traces -> %d elution peaks -> %d features", len(mass_traces), len(elution_peaks), features.size())
    return features


def find_features(
    mzml: str | Path,
    features_out: str | Path,
    precursors_out: str | Path,
    instrument: InstrumentSettings,
    settings: FeatureFindingSettings,
) -> oms.FeatureMap:
    """
    Detect the features of one run and correct its MS2 precursors.

    Writes the features (featureXML) and the corrected precursor m/z of every
    MS2 spectrum (TSV, see ``msdata``). The mzML file itself is not modified.
    """
    experiment = load_experiment(mzml)
    check_centroided(experiment, mzml)

    # Step 1: precursor m/z -> most intense MS1 peak within the tolerance.
    # (The three lists are optional diagnostic outputs, unused here.)
    oms.PrecursorCorrection.correctToHighestIntensityMS1Peak(
        experiment, settings.precursor_peak_tolerance_ppm, True, [], [], []
    )

    # Step 2: feature detection.
    features = detect_features(experiment, instrument, settings)
    features.setPrimaryMSRunPath([str(mzml).encode()])

    # Step 3: precursor m/z -> m/z of the feature that contains it. Arguments:
    # rt tolerance 0 s (the MS2 must fall inside the feature), m/z tolerance in
    # ppm, don't trust the reported charge, don't keep the original precursor,
    # only the nearest feature, and look at the first 2 isotope traces.
    corrected = oms.PrecursorCorrection.correctToNearestFeature(
        features, experiment, 0.0, settings.precursor_feature_tolerance_ppm, True, False, False, False, 2, 0
    )
    precursors = precursor_table(experiment)
    log.info("%d MS2 spectra, %d precursors corrected to a feature", len(precursors), len(corrected))

    store_feature_map(features_out, features)
    precursors.to_csv(precursors_out, sep="\t", index=False)
    return features
