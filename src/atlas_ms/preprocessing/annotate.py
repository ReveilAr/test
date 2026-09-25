"""
Per-run feature annotation: adduct grouping and MS2 mapping.

1. Adduct grouping (UmetaFlow's `adduct_annotations_FFM` rule): one molecule
   often gives several ions at the same RT, e.g. [M+H]+ and [M+Na]+.
   MetaboliteAdductDecharger finds features whose m/z differences match the
   configured adducts and groups them. These groups become the "ion identity"
   links of Ion Identity Molecular Networking (IIMN).
2. MS2 mapping (UmetaFlow's `IDMapper_FFM` rule): each MS2 spectrum is
   attached to the feature whose convex hull contains its precursor m/z and
   RT. Downstream, the MGF export takes each feature's spectra from these
   attachments. OpenMS stores them as (empty) peptide identifications, which
   is why this step uses IDMapper.
"""

import logging
from pathlib import Path

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


def attach_ms2(features: oms.FeatureMap, experiment: oms.MSExperiment, run_name: str) -> None:
    """
    Attach every MS2 spectrum of ``experiment`` to the feature that contains
    its precursor (in place). Both must be on the same (aligned) RT axis.
    """
    # IDMapper records the attachments under an identification "run". It must
    # have a unique name per sample, or the consensus map can't be saved later.
    id_run = oms.ProteinIdentification()
    id_run.setIdentifier(run_name)
    # use_centroid_rt=False: match against the feature's RT extent (convex
    # hull), not only its apex; use_centroid_mz=True: match on the feature
    # m/z. These are the IDMapper tool defaults used by UmetaFlow.
    oms.IDMapper().annotate(features, oms.PeptideIdentificationList(), [id_run], False, True, experiment)
    drop_empty_ms2(features, experiment)


def drop_empty_ms2(features: oms.FeatureMap, experiment: oms.MSExperiment) -> None:
    """
    Detach MS2 spectra that have no peaks (in place).

    Instruments record some MS2 scans without any centroid, e.g. scans
    triggered on noise. They are useless, and harmful at export: the GNPS MGF
    writer takes each feature's spectrum from the run where the feature is
    most intense, and skips the whole feature if that spectrum is empty,
    even when other runs have good spectra. A feature whose only spectra are
    empty becomes a feature without MS2.
    """
    n_dropped = 0
    for index in range(features.size()):
        feature = features[index]  # a copy (pyOpenMS)
        attached = feature.getPeptideIdentifications()
        usable = oms.PeptideIdentificationList()
        for ms2 in attached:
            if experiment[int(ms2.getMetaValue("spectrum_index"))].size() > 0:
                usable.push_back(ms2)
        if usable.size() < attached.size():
            n_dropped += attached.size() - usable.size()
            feature.setPeptideIdentifications(usable)
            features[index] = feature  # write the modified copy back
    if n_dropped:
        log.info("%d MS2 spectra without peaks were not attached to features", n_dropped)


def annotate_features(
    features: oms.FeatureMap,
    experiment: oms.MSExperiment,
    run_name: str,
    adducts: AdductSettings,
) -> oms.FeatureMap:
    """Adduct grouping followed by MS2 mapping; returns the annotated features."""
    grouped = group_adducts(features, adducts)
    attach_ms2(grouped, experiment, run_name)
    # Keep the reference to the run's mzML file (pyOpenMS fills a list in place).
    run_paths = []
    features.getPrimaryMSRunPath(run_paths)
    grouped.setPrimaryMSRunPath(run_paths)
    with_ms2 = sum(1 for feature in grouped if feature.getPeptideIdentifications().size() > 0)
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
