"""
Gap filling: targeted re-extraction of features missing in some runs
(UmetaFlow's requantification rules).

A compound can be missed by untargeted feature finding in a run where it is
weak, which leaves gaps in the feature table. UmetaFlow's scheme, kept here:

1. ``plan_gap_filling``: split the linked consensus features into
   * "complete" ones, found in every run: their features are kept as they are;
   * "incomplete" ones: they become targets (m/z, charge, RT) for
     re-extraction.
2. ``fill_gaps`` (per run): FeatureFinderMetaboIdent extracts every target
   from the raw data, in *every* run, so all values of an incomplete feature
   come from the same method. The extracted features are merged with the
   complete ones, then adduct grouping and MS2 mapping are redone.
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


def plan_gap_filling(consensus_file: str | Path, targets_out: str | Path, complete_out: str | Path) -> None:
    """
    Write the re-extraction targets and the list of features to keep.

    targets.tsv:    target, mz, charge, rt       (one row per incomplete consensus feature)
    complete.tsv:   map_index, feature_id        (features of the complete consensus features)
    """
    consensus = load_consensus_map(consensus_file)
    n_runs = len(consensus.getColumnHeaders())
    targets, complete = [], []
    for index, cfeature in enumerate(consensus):
        if cfeature.size() == n_runs:
            # Unique ids are 64-bit unsigned integers: stored as text to avoid overflow.
            complete.extend((h.getMapIndex(), str(h.getUniqueId())) for h in cfeature.getFeatureList())
        else:
            targets.append((f"target_{index}", cfeature.getMZ(), cfeature.getCharge(), cfeature.getRT()))
    pd.DataFrame(targets, columns=["target", "mz", "charge", "rt"]).to_csv(targets_out, sep="\t", index=False)
    pd.DataFrame(complete, columns=["map_index", "feature_id"]).to_csv(complete_out, sep="\t", index=False)
    log.info("%d complete consensus features, %d targets for gap filling", consensus.size() - len(targets), len(targets))


def extract_targets(
    experiment: oms.MSExperiment,
    targets: pd.DataFrame,
    polarity: str,
    instrument: InstrumentSettings,
    settings: GapFillingSettings,
    mzml: str | Path,
) -> oms.FeatureMap:
    """Targeted extraction of ``targets`` from one run (FeatureFinderMetaboIdent)."""
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


def fill_gaps(
    mzml: str | Path,
    precursors: str | Path,
    features_in: str | Path,
    trafo: str | Path,
    targets_file: str | Path,
    complete_file: str | Path,
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

    # Keep the features that belong to complete consensus features.
    complete = pd.read_csv(complete_file, sep="\t", dtype={"feature_id": str})
    keep_ids = set(complete.loc[complete["map_index"] == map_index, "feature_id"])
    merged = oms.FeatureMap(features)
    merged.clear(False)  # False: empty the features but keep the map's metadata
    for feature in features:
        if str(feature.getUniqueId()) in keep_ids:
            merged.push_back(feature)
    n_kept = merged.size()

    # Re-extract every target in this run and add the results.
    targets = pd.read_csv(targets_file, sep="\t")
    for feature in extract_targets(experiment, targets, adducts.polarity, instrument, settings, mzml):
        merged.push_back(feature)
    merged.setUniqueIds()
    log.info("%s: %d features kept + %d re-extracted (%d targets)", run_name, n_kept, merged.size() - n_kept, len(targets))

    store_feature_map(features_out, annotate_features(merged, experiment, run_name, adducts))
