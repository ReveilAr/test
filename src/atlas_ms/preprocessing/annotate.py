"""
Per-run feature annotation: adduct grouping and MS2 mapping.

1. Adduct grouping (UmetaFlow's `adduct_annotations_FFM` rule): one molecule
   often gives several ions at the same RT, e.g. [M+H]+ and [M+Na]+.
   MetaboliteAdductDecharger finds features whose m/z differences match the
   configured adducts and groups them. These groups become the "ion identity"
   links of Ion Identity Molecular Networking (IIMN).
2. MS2 mapping (UmetaFlow's `IDMapper_FFM` rule): each MS2 spectrum is
   attached to the feature it was recorded on, i.e. the feature whose
   retention-time extent contains the spectrum and whose m/z is the
   precursor m/z. The spectrum indices are stored on the feature (meta
   value "ms2_spectra"); the MGF export takes each feature's spectrum from
   them.

   UmetaFlow does this with OpenMS's IDMapper, a proteomics tool that stores
   the attachments as empty *peptide* identifications. ``attach_ms2`` is the
   same matching rule (same tolerances), written for small molecules, with
   three differences:

   * a spectrum whose precursor fits several features goes to the closest
     one in m/z (IDMapper: the first one in the map);
   * an unknown charge (0, frequent for the precursors of small molecules)
     matches any charge (IDMapper requires equal charges and loses these
     spectra);
   * MS2 spectra without any peak are not attached (see ``attach_ms2``).
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pyopenms as oms

from atlas_ms.config import AdductSettings
from atlas_ms.preprocessing.msdata import (
    apply_trafo,
    load_feature_map,
    load_run,
    load_trafo,
    set_parameters,
    store_feature_map,
)

log = logging.getLogger(__name__)


def group_adducts(features: oms.FeatureMap, adducts: AdductSettings) -> oms.FeatureMap:
    """
    Group the features that are different adducts of the same molecule.
    Returns a new FeatureMap. Each feature is annotated with its most likely
    adduct and a group id ("dc_charge_adducts", "Group" meta values).
    """
    if adducts.polarity == "positive":
        # Singly charged positive ions only, at most one neutral loss/gain.
        charges = {"charge_min": 1, "charge_max": 1, "charge_span_max": 1}
    else:
        # UmetaFlow's negative-mode settings.
        charges = {"charge_min": -2, "charge_max": 0, "charge_span_max": 3, "negative_mode": "true"}

    decharger = oms.MetaboliteFeatureDeconvolution()
    set_parameters(decharger, {
        "potential_adducts": adducts.openms_adducts(),
        "max_neutrals": 1,
        "retention_max_diff": 3.0,  # s: adducts of one molecule co-elute
        "retention_max_diff_local": 3.0,
        **charges,
    })
    grouped = oms.FeatureMap()
    decharger.compute(features, grouped, oms.ConsensusMap(), oms.ConsensusMap())
    return grouped


# Matching tolerances of an MS2 precursor to a feature (IDMapper's defaults,
# used by UmetaFlow).
MS2_RT_TOLERANCE_S = 5.0  # added on both sides of the feature's RT extent
MS2_MZ_TOLERANCE_PPM = 20.0


def feature_extents(features: oms.FeatureMap) -> pd.DataFrame:
    """
    Retention-time extent (start, end), m/z and charge of every feature. The
    extent is that of the feature's convex hull (all its mass traces), or
    its apex when it has no hull.
    """
    rows = []
    for feature in features:
        if feature.getConvexHulls():
            box = feature.getConvexHull().getBoundingBox()
            start, end = box.minPosition()[0], box.maxPosition()[0]
        else:
            start = end = feature.getRT()
        rows.append((start, end, feature.getMZ(), abs(feature.getCharge())))
    return pd.DataFrame(rows, columns=["rt_start", "rt_end", "mz", "charge"])


def attach_ms2(features: oms.FeatureMap, experiment: oms.MSExperiment) -> None:
    """
    Attach every MS2 spectrum of ``experiment`` to the feature it was
    recorded on (in place). Both must be on the same (aligned) RT axis.

    Candidates are the features whose RT extent, widened by
    ``MS2_RT_TOLERANCE_S`` on each side, contains the spectrum, whose m/z is
    within ``MS2_MZ_TOLERANCE_PPM`` of the precursor m/z, and whose charge
    agrees (0 = unknown agrees with anything). The spectrum goes to the
    candidate closest in m/z. The spectrum indices (positions in the run)
    are stored in the feature's "ms2_spectra" meta value.

    MS2 spectra without peaks are not attached. Instruments record some,
    e.g. scans triggered on noise. They are useless, and a feature whose
    only spectra are empty is a feature without MS2.
    """
    extents = feature_extents(features)
    start = extents["rt_start"].to_numpy() - MS2_RT_TOLERANCE_S
    end = extents["rt_end"].to_numpy() + MS2_RT_TOLERANCE_S
    mz = extents["mz"].to_numpy()
    charge = extents["charge"].to_numpy()

    attached = {}  # feature position -> spectrum indices
    n_ms2 = n_empty = n_unassigned = 0
    for index, spectrum in enumerate(experiment):
        if spectrum.getMSLevel() != 2 or not spectrum.getPrecursors():
            continue
        n_ms2 += 1
        if spectrum.size() == 0:
            n_empty += 1
            continue
        precursor = spectrum.getPrecursors()[0]
        rt, precursor_charge = spectrum.getRT(), abs(precursor.getCharge())
        ppm = np.abs(mz - precursor.getMZ()) / mz * 1e6
        candidates = ((start <= rt) & (rt <= end) & (ppm <= MS2_MZ_TOLERANCE_PPM)
                      & ((charge == precursor_charge) | (charge == 0) | (precursor_charge == 0)))
        if not candidates.any():
            n_unassigned += 1
            continue
        closest = np.flatnonzero(candidates)[np.argmin(ppm[candidates])]
        attached.setdefault(int(closest), []).append(index)

    for position, indices in attached.items():
        feature = features[position]  # a copy (pyOpenMS)
        feature.setMetaValue("ms2_spectra", indices)
        features[position] = feature  # write the modified copy back
    log.info("%d MS2 spectra: %d attached to %d features, %d outside any feature, %d without peaks (skipped)",
             n_ms2, n_ms2 - n_empty - n_unassigned, len(attached), n_unassigned, n_empty)


def ms2_spectra(feature) -> list[int]:
    """The run's spectrum indices of the MS2 spectra attached to a feature (``attach_ms2``)."""
    return list(feature.getMetaValue("ms2_spectra")) if feature.metaValueExists("ms2_spectra") else []


def annotate_features(
    features: oms.FeatureMap,
    experiment: oms.MSExperiment,
    run_name: str,
    adducts: AdductSettings,
) -> oms.FeatureMap:
    """Adduct grouping followed by MS2 mapping; returns the annotated features."""
    grouped = group_adducts(features, adducts)
    attach_ms2(grouped, experiment)
    # Keep the reference to the run's mzML file (pyOpenMS fills a list in place).
    run_paths = []
    features.getPrimaryMSRunPath(run_paths)
    grouped.setPrimaryMSRunPath(run_paths)
    with_ms2 = sum(1 for feature in grouped if feature.metaValueExists("ms2_spectra"))
    log.info("%s: %d features, %d with MS2", run_name, grouped.size(), with_ms2)
    return grouped


def annotate_run(
    mzml: str | Path,
    precursors: str | Path,
    features_in: str | Path,
    trafo: str | Path,
    features_out: str | Path,
    run_name: str,
    adducts: AdductSettings,
) -> None:
    """File-level entry point: align one run's features and spectra, then annotate."""
    experiment = load_run(mzml, precursors, trafo)
    features = load_feature_map(features_in)
    apply_trafo(features, load_trafo(trafo))
    store_feature_map(features_out, annotate_features(features, experiment, run_name, adducts))
