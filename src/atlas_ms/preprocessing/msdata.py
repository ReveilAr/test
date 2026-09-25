"""
Reading and writing LC-MS data with pyOpenMS, and the small "correction"
files that describe how a run must be adjusted.

The pipeline never writes a modified copy of an mzML file. On a laptop,
100 runs would mean tens of GB of duplicates. Two small files per run
record the corrections instead:

* ``work/features/<sample>.precursors.tsv``: the corrected precursor m/z of
  every MS2 spectrum (written by feature finding);
* ``work/alignment/<sample>.trafoXML``: the retention-time transformation
  found by the alignment.

Every step that needs the spectra loads the original mzML and re-applies
both in memory with ``load_run`` (fast C++ operations).
"""

from pathlib import Path

import pandas as pd
import pyopenms as oms


# ---- mzML ----

def load_experiment(path: str | Path) -> oms.MSExperiment:
    """Load a whole mzML file (all MS levels) into memory."""
    experiment = oms.MSExperiment()
    oms.MzMLFile().load(str(path), experiment)
    return experiment


def check_centroided(experiment: oms.MSExperiment, path: str | Path) -> None:
    """
    Raise if the MS1 data is in profile mode. Feature finding needs centroids,
    and ThermoRawFileParser centroids by default. An mzML made by another
    tool must be converted with peak picking.
    """
    for spectrum in experiment:
        if spectrum.getMSLevel() == 1:
            if spectrum.getType() == oms.SpectrumSettings.SpectrumType.PROFILE:
                raise ValueError(f"{path} contains profile MS1 spectra; convert it with centroiding (peak picking)")
            return


def ms1_experiment(experiment: oms.MSExperiment, min_intensity: float = 0.0) -> oms.MSExperiment:
    """
    Copy of ``experiment`` with only its MS1 scans.

    With ``min_intensity`` > 0, scans without any peak at or above that
    intensity are skipped as well. This works around an OpenMS 3.5 issue:
    MassTraceDetection loses most traces when some scans have no peak above
    its noise threshold. Skipping those scans cannot remove a trace point,
    because any such point would be a peak above the threshold.
    """
    keep = []
    for spectrum in experiment:
        if spectrum.getMSLevel() != 1:
            continue
        if min_intensity > 0:
            _, intensities = spectrum.get_peaks()
            if intensities.size == 0 or intensities.max() < min_intensity:
                continue
        keep.append(spectrum)
    ms1 = oms.MSExperiment()
    ms1.setSpectra(keep)
    ms1.updateRanges()
    return ms1


# ---- precursor corrections ----

def precursor_table(experiment: oms.MSExperiment) -> pd.DataFrame:
    """
    Current precursor m/z of every MS2 spectrum. ``index`` is the spectrum's
    position in the file; ``native_id`` is kept to check that a table is
    applied to the file it was made from.
    """
    rows = [
        (index, spectrum.getNativeID(), spectrum.getPrecursors()[0].getMZ())
        for index, spectrum in enumerate(experiment)
        if spectrum.getMSLevel() == 2 and spectrum.getPrecursors()
    ]
    return pd.DataFrame(rows, columns=["index", "native_id", "precursor_mz"])


def apply_precursor_table(experiment: oms.MSExperiment, table: pd.DataFrame) -> None:
    """Set the precursor m/z of the MS2 spectra listed in ``table`` (in place)."""
    for index, native_id, mz in table.itertuples(index=False):
        spectrum = experiment[index]  # pyOpenMS returns a copy...
        if spectrum.getNativeID() != native_id:
            raise ValueError(f"Precursor table does not match the mzML file (spectrum {index} is not {native_id})")
        precursors = spectrum.getPrecursors()
        precursors[0].setMZ(mz)
        spectrum.setPrecursors(precursors)
        experiment[index] = spectrum  # ...so the modified copy is written back


# ---- retention-time transformations ----

def load_trafo(path: str | Path) -> oms.TransformationDescription:
    trafo = oms.TransformationDescription()
    oms.TransformationXMLFile().load(str(path), trafo, True)  # True: fit the stored model
    return trafo


def store_trafo(path: str | Path, trafo: oms.TransformationDescription) -> None:
    oms.TransformationXMLFile().store(str(path), trafo)


def apply_trafo(data, trafo: oms.TransformationDescription) -> None:
    """
    Move the retention times of an MSExperiment or FeatureMap (in place) onto
    the common, aligned time axis. The original times are kept as the
    "original_RT" meta value.
    """
    oms.MapAlignmentTransformer().transformRetentionTimes(data, trafo, True)


# ---- feature and consensus maps ----

def load_feature_map(path: str | Path) -> oms.FeatureMap:
    features = oms.FeatureMap()
    oms.FeatureXMLFile().load(str(path), features)
    return features


def store_feature_map(path: str | Path, features: oms.FeatureMap) -> None:
    oms.FeatureXMLFile().store(str(path), features)


def load_consensus_map(path: str | Path) -> oms.ConsensusMap:
    consensus = oms.ConsensusMap()
    oms.ConsensusXMLFile().load(str(path), consensus)
    return consensus


def store_consensus_map(path: str | Path, consensus: oms.ConsensusMap) -> None:
    oms.ConsensusXMLFile().store(str(path), consensus)


# ---- one run, ready to use ----

def load_run(mzml: str | Path, precursors: str | Path, trafo: str | Path) -> oms.MSExperiment:
    """
    Load a run with both corrections applied: precursor m/z from feature
    finding, and retention times on the aligned axis.
    """
    experiment = load_experiment(mzml)
    apply_precursor_table(experiment, pd.read_csv(precursors, sep="\t"))
    apply_trafo(experiment, load_trafo(trafo))
    return experiment


def set_parameters(algorithm, values: dict) -> None:
    """Set some parameters of a pyOpenMS algorithm, keeping defaults for the others."""
    parameters = algorithm.getDefaults()
    for key, value in values.items():
        parameters.setValue(key, value)
    algorithm.setParameters(parameters)
