"""
Gap filling: targeted re-extraction of features missing in some runs.

A compound can be missed by untargeted feature finding in a run where it is
weak, which leaves gaps in the feature table. Gap filling looks for it again,
in a targeted way, where it is missing:

1. ``plan_gap_filling``: every consensus feature not found in all runs becomes
   a target (m/z, charge, RT and its measured isotope pattern).
2. ``fill_gaps`` (per run): FeatureFinderMetaboIdent extracts every target
   from the raw data. Then:
   * features found by untargeted feature finding are always kept as they
     are: re-extraction never replaces a detected value (it could pick a
     neighbouring isomer, which is common among lipids);
   * a re-extracted feature is added only where the run had no feature for
     that target: that is, it fills a gap;
   * re-extracted intensities are put on the scale of the detected ones (see
     ``intensity_scale``), so that every value of the table means the same
     thing. This matters for lipid class sums as much as for single features.
   Adduct grouping and MS2 mapping are then redone on the merged features.
3. The runs are linked again (``linking.link_runs``).

This differs from UmetaFlow, which discards the detected values of
incomplete features and re-extracts them in every run. That loses values
when re-extraction fails, and mixes two intensity definitions between
features.
"""

import logging
from pathlib import Path

import numpy as np
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

# Fewer detected/re-extracted pairs than this in a run: the scale factor
# would be unreliable, so re-extracted intensities are left unscaled.
MIN_PAIRS_FOR_SCALE = 10


# --------------------------------------------------------------------------
# Planning (all runs)
# --------------------------------------------------------------------------

def isotope_patterns(feature_files: list[str | Path]) -> dict[tuple[int, str], list[float]]:
    """
    Measured isotope pattern of every feature: (map index, feature id) ->
    relative intensities of its mass traces (M, M+1, M+2...), summing to 1.
    FeatureFinderMetabo stores them as the "masstrace_intensity" meta value.
    """
    patterns = {}
    for map_index, path in enumerate(feature_files):
        for feature in load_feature_map(path):
            traces = np.asarray(feature.getMetaValue("masstrace_intensity"), dtype=float)
            if traces.size and traces.sum() > 0:
                patterns[(map_index, str(feature.getUniqueId()))] = list(traces / traces.sum())
    return patterns


def plan_gap_filling(
    consensus_file: str | Path,
    feature_files: list[str | Path],
    targets_out: str | Path,
    members_out: str | Path,
) -> None:
    """
    Write the re-extraction targets, and where every existing feature belongs.

    targets.tsv:  target, mz, charge, rt, isotopes   one row per incomplete consensus
                  (isotopes: measured pattern,        feature
                  ";"-separated)
    members.tsv:  map_index, feature_id, target      one row per feature; target is
                                                     empty for complete consensus features

    ``feature_files`` are the feature maps that were linked (in map-index
    order); they provide the isotope patterns.
    """
    consensus = load_consensus_map(consensus_file)
    patterns = isotope_patterns(feature_files)
    n_runs = len(consensus.getColumnHeaders())
    targets, members = [], []
    for index, cfeature in enumerate(consensus):
        handles = cfeature.getFeatureList()
        target = ""
        if cfeature.size() < n_runs:
            target = f"target_{index}"
            # Isotope pattern of the most intense occurrence: the best measured one.
            best = max(handles, key=lambda h: h.getIntensity())
            pattern = patterns.get((best.getMapIndex(), str(best.getUniqueId())), [])
            isotopes = ";".join(f"{p:.4f}" for p in pattern)
            targets.append((target, cfeature.getMZ(), cfeature.getCharge(), cfeature.getRT(), isotopes))
        # Unique ids are 64-bit unsigned integers: stored as text to avoid overflow.
        members.extend((h.getMapIndex(), str(h.getUniqueId()), target) for h in handles)
    pd.DataFrame(targets, columns=["target", "mz", "charge", "rt", "isotopes"]).to_csv(targets_out, sep="\t", index=False)
    pd.DataFrame(members, columns=["map_index", "feature_id", "target"]).to_csv(members_out, sep="\t", index=False)
    log.info("%d complete consensus features, %d targets for gap filling", consensus.size() - len(targets), len(targets))


# --------------------------------------------------------------------------
# Re-extraction (one run)
# --------------------------------------------------------------------------

def read_members(members_file: str | Path, map_index: int) -> dict[str, str]:
    """Feature id -> target name ("" for complete features), for one run."""
    members = pd.read_csv(members_file, sep="\t", dtype=str, keep_default_na=False)
    members = members[members["map_index"] == str(map_index)]
    return dict(zip(members["feature_id"], members["target"]))


def monoisotopic_share(feature: oms.Feature) -> float:
    """
    Fraction of a re-extracted feature's intensity that belongs to its
    monoisotopic trace. FeatureFinderMetaboIdent sums all extracted isotope
    traces, while FeatureFinderMetabo reports the monoisotopic one only.
    """
    for trace in feature.getSubordinates():
        if not trace.metaValueExists("native_id"):
            continue
        native_id = trace.getMetaValue("native_id")
        native_id = native_id.decode() if isinstance(native_id, bytes) else str(native_id)
        if native_id.endswith("_i0"):  # trace of isotope 0
            return float(trace.getMetaValue("isotope_probability"))
    return 1.0  # no isotope information: keep the intensity as it is


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
    Each extracted feature carries its target name as the "label" meta value,
    and its intensity is that of its monoisotopic trace.
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
        # The measured isotope pattern. Without a formula, OpenMS would
        # otherwise fall back to a peptide model (and log an error per target).
        isotopes = [float(p) for p in str(target.isotopes).split(";") if p not in ("", "nan")] or [0.0]
        # Arguments: name, formula (unknown), mass, charges, RTs, RT ranges
        # (0 = use the extraction window), isotope pattern.
        library.append(oms.FeatureFinderMetaboIdentCompound(
            target.target, "", mass, [charge], [target.rt], [0.0], isotopes
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

    # Keep only the monoisotopic part of each intensity, as feature finding does.
    for index in range(features.size()):
        feature = features[index]  # a copy (pyOpenMS)
        feature.setIntensity(feature.getIntensity() * monoisotopic_share(feature))
        features[index] = feature
    return features


def intensity_scale(originals: dict[str, float], extracted: dict[str, float]) -> tuple[float, int]:
    """
    Factor that converts re-extracted intensities to the scale of detected ones.

    Both are monoisotopic peak areas, but computed differently (trace area
    vs. fitted model), so they differ by a roughly constant factor. It is
    measured on the targets that were both detected and re-extracted in the
    run, as the median ratio, which is robust to the few pairs where
    re-extraction picked another peak. ``originals`` and ``extracted`` map
    target -> intensity. Returns (factor, number of pairs).
    """
    common = [t for t in sorted(set(originals) & set(extracted)) if extracted[t] > 0]
    if len(common) < MIN_PAIRS_FOR_SCALE:
        return 1.0, len(common)
    detected = np.array([originals[t] for t in common])
    ratios = detected / np.array([extracted[t] for t in common])

    # Diagnostic: the gaps are mostly weak features, so the factor is only
    # right for them if the ratio does not depend on intensity. Compare the
    # median ratio of the weakest and strongest quarters of the pairs.
    q1, q3 = np.percentile(detected, [25, 75])
    log.info(
        "Scale from %d pairs: median %.3f (IQR %.3f-%.3f); weakest quarter %.3f, strongest quarter %.3f",
        len(ratios), np.median(ratios), *np.percentile(ratios, [25, 75]),
        np.median(ratios[detected <= q1]), np.median(ratios[detected >= q3]),
    )
    return float(np.median(ratios)), len(ratios)


def merge_gap_filled(
    features: oms.FeatureMap,
    members: dict[str, str],
    extracted: oms.FeatureMap,
) -> tuple[oms.FeatureMap, int, float]:
    """
    The features of one run after gap filling:

    * every feature found by feature finding that belongs to a consensus
      feature (``members``), unchanged;
    * re-extracted features of the targets this run had no feature for,
      with their intensity rescaled (``intensity_scale``).

    Returns (merged features, number of gaps filled, scale factor used).
    """
    # Target -> intensity, for the detected and the re-extracted features.
    originals = {}
    for feature in features:
        target = members.get(str(feature.getUniqueId()))
        if target:
            originals[target] = feature.getIntensity()
    extracted_intensity = {f.getMetaValue("label"): f.getIntensity() for f in extracted}
    scale, n_pairs = intensity_scale(originals, extracted_intensity)
    if n_pairs < MIN_PAIRS_FOR_SCALE:
        log.warning("Only %d detected/re-extracted pairs: re-extracted intensities left unscaled", n_pairs)

    merged = oms.FeatureMap(features)
    merged.clear(False)  # False: empty the features but keep the map's metadata
    for feature in features:
        if str(feature.getUniqueId()) in members:
            merged.push_back(feature)
    n_filled = 0
    for feature in extracted:
        if feature.getMetaValue("label") not in originals:
            feature.setIntensity(feature.getIntensity() * scale)
            # Mark it: the export flags gap-filled values (they are less
            # precise than detected ones, see intensity_scale).
            feature.setMetaValue("gap_filled", "true")
            merged.push_back(feature)
            n_filled += 1
    merged.setUniqueIds()
    return merged, n_filled, scale


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

    targets = pd.read_csv(targets_file, sep="\t", keep_default_na=False)
    extracted = extract_targets(experiment, targets, adducts.polarity, instrument, settings, mzml)
    merged, n_filled, scale = merge_gap_filled(features, members, extracted)
    log.info(
        "%s: %d targets, %d re-extracted, %d gaps filled (intensity scale %.3f), %d detected features kept",
        run_name, len(targets), extracted.size(), n_filled, scale, merged.size() - n_filled,
    )
    store_feature_map(features_out, annotate_features(merged, experiment, run_name, adducts))
