"""
Gap filling: targeted re-extraction of features missing in some runs
(UmetaFlow's requantification rules, plus one safeguard).

A compound can be missed by untargeted feature finding in a run where it is
weak, which leaves gaps in the feature table. The scheme:

1. ``plan_gap_filling``: split the linked consensus features into
   * "complete" ones, found in every run: their features are kept as they are;
   * "incomplete" ones: they become targets (m/z, charge, RT) for
     re-extraction.
2. ``fill_gaps`` (per run): FeatureFinderMetaboIdent extracts every target
   from the raw data, in *every* run, so that all values of an incomplete
   feature come from the same method (UmetaFlow's choice).
   Safeguard (not in UmetaFlow): when the extraction of a target fails in a
   run where untargeted feature finding had found it, the original feature
   is kept. Otherwise that value is lost, and a feature whose extraction
   fails everywhere disappears from the table.
   Adduct grouping and MS2 mapping are then redone on the merged features.
3. The runs are linked again (``linking.link_runs``).
"""

import logging
from pathlib import Path

import pandas as pd
import pyopenms as oms

from atlas_ms.config import AdductSettings, GapFillingSettings, InstrumentSettings
from atlas_ms.preprocessing.annotate import annotate_features
from atlas_ms.preprocessing.msdata import (
    apply_trafo,
    load_consensus_map,
    load_feature_map,
    load_run,
    load_trafo,
    ms1_experiment,
    set_parameters,
    store_feature_map,
)

log = logging.getLogger(__name__)

PROTON_MASS = 1.007276466  # u


def plan_gap_filling(consensus_file: str | Path, targets_out: str | Path, members_out: str | Path) -> None:
    """
    Write the re-extraction targets, and where every existing feature belongs.

    targets.tsv:  target, mz, charge, rt           one row per incomplete consensus feature
    members.tsv:  map_index, feature_id, target    one row per feature of any consensus
                                                   feature; target is empty for complete ones
    """
    consensus = load_consensus_map(consensus_file)
    n_runs = len(consensus.getColumnHeaders())
    targets, members = [], []
    for index, cfeature in enumerate(consensus):
        target = ""
        if cfeature.size() < n_runs:
            target = f"target_{index}"
            targets.append((target, cfeature.getMZ(), cfeature.getCharge(), cfeature.getRT()))
        # Unique ids are 64-bit unsigned integers: stored as text to avoid overflow.
        members.extend((h.getMapIndex(), str(h.getUniqueId()), target) for h in cfeature.getFeatureList())
    pd.DataFrame(targets, columns=["target", "mz", "charge", "rt"]).to_csv(targets_out, sep="\t", index=False)
    pd.DataFrame(members, columns=["map_index", "feature_id", "target"]).to_csv(members_out, sep="\t", index=False)
    log.info("%d complete consensus features, %d targets for gap filling", consensus.size() - len(targets), len(targets))


def read_members(members_file: str | Path, map_index: int) -> dict[str, str]:
    """Feature id -> target name ("" for complete features), for one run."""
    members = pd.read_csv(members_file, sep="\t", dtype=str, keep_default_na=False)
    members = members[members["map_index"] == str(map_index)]
    return dict(zip(members["feature_id"], members["target"]))


def extract_targets(
    experiment: oms.MSExperiment,
    targets: pd.DataFrame,
    polarity: str,
    instrument: InstrumentSettings,
    settings: GapFillingSettings,
    mzml: str | Path,
) -> oms.FeatureMap:
    """
    Targeted extraction of ``targets`` from one run (FeatureFinderMetaboIdent).
    Each extracted feature carries its target name as the "label" meta value.
    """
    features = oms.FeatureMap()
    if targets.empty:
        return features

    sign = 1 if polarity == "positive" else -1
    library = []
    for target in targets.itertuples():
        # The linker reports charge 0 when unknown: assume a singly charged ion.
        charge = sign * max(abs(int(target.charge)), 1)
        # FeatureFinderMetaboIdent expects a neutral mass and recomputes the
        # m/z as (mass + charge * proton) / |charge|. Any mass that gives back
        # the observed m/z works, even if the ion is not a protonated molecule.
        mass = target.mz * abs(charge) - charge * PROTON_MASS
        # Arguments: name, formula (unknown), mass, charges, RTs, RT ranges
        # (0 = use the extraction window), isotope distribution (0 = compute).
        library.append(oms.FeatureFinderMetaboIdentCompound(
            target.target, "", mass, [charge], [target.rt], [0.0], [0.0]
        ))

    finder = oms.FeatureFinderAlgorithmMetaboIdent()
    finder.setMSData(ms1_experiment(experiment))
    set_parameters(finder, {
        "extract:mz_window": settings.mz_window_ppm,
        "extract:rt_window": settings.rt_window_s,
        # Width of the smoothing applied before peak detection: the expected
        # peak width. (The OpenMS default of 60 s suits long gradients only.)
        "detect:peak_width": instrument.chrom_fwhm_s,
    })
    finder.run(library, features, str(mzml))
    return features


def merge_gap_filled(
    features: oms.FeatureMap,
    members: dict[str, str],
    extracted: oms.FeatureMap,
) -> oms.FeatureMap:
    """
    The features of one run after gap filling:

    * features of complete consensus features (``members`` target "");
    * every re-extracted feature;
    * original features of incomplete consensus features whose target was
      *not* re-extracted in this run (the safeguard, see module docstring).

    Features that belong to no consensus feature are dropped. The linker
    places every feature in a consensus feature, so none are expected.
    """
    extracted_targets = {feature.getMetaValue("label") for feature in extracted}
    merged = oms.FeatureMap(features)
    merged.clear(False)  # False: empty the features but keep the map's metadata
    for feature in features:
        target = members.get(str(feature.getUniqueId()))
        if target == "" or (target and target not in extracted_targets):
            merged.push_back(feature)
    for feature in extracted:
        merged.push_back(feature)
    merged.setUniqueIds()
    return merged


def fill_gaps(
    mzml: str | Path,
    precursors: str | Path,
    features_in: str | Path,
    trafo: str | Path,
    targets_file: str | Path,
    members_file: str | Path,
    map_index: int,
    features_out: str | Path,
    run_name: str,
    instrument: InstrumentSettings,
    adducts: AdductSettings,
    settings: GapFillingSettings,
) -> None:
    """
    Gap filling for one run.

    ``features_in`` are the run's features straight out of feature finding,
    before adduct grouping, because grouping and MS2 mapping are redone on
    the merged set. ``map_index`` is the position of the run in the sample
    table.
    """
    experiment = load_run(mzml, precursors, trafo)
    features = load_feature_map(features_in)
    apply_trafo(features, load_trafo(trafo))
    members = read_members(members_file, map_index)

    targets = pd.read_csv(targets_file, sep="\t")
    extracted = extract_targets(experiment, targets, adducts.polarity, instrument, settings, mzml)
    merged = merge_gap_filled(features, members, extracted)

    n_complete = sum(1 for target in members.values() if target == "")
    n_rescued = merged.size() - n_complete - extracted.size()
    log.info(
        "%s: %d targets, %d re-extracted, %d original features kept where re-extraction failed, "
        "%d features from complete consensus features",
        run_name, len(targets), extracted.size(), n_rescued, n_complete,
    )
    store_feature_map(features_out, annotate_features(merged, experiment, run_name, adducts))
