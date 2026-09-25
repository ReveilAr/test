"""
Final feature table and GNPS/FBMN export (UmetaFlow's `GNPS_export` rules).

Outputs:

* ``results/features.parquet``: one row per feature (m/z, RT, MS2, adducts...)
* ``results/quant.parquet``: feature intensities, one column per sample
* ``results/quant_gap_filled.parquet``: same shape, True where the value was
  re-extracted by gap filling rather than detected
* ``results/gnps/``: the four files of GNPS Feature-Based Molecular
  Networking with the "OpenMS" input format (MS2 spectra as MGF,
  quantification table, sample metadata, IIMN adduct pairs). They are what
  downstream tools (FBMN-STATS-style statistics, GNPS2) expect, so the
  format is not reinvented.

Feature ids: GNPS numbers features 1..N in the order of the consensus map
("row ID" in the quantification table, "SCANS" in the MGF). Only features
with MS2 go to GNPS, but the internal tables keep all features, because
features without MS2 still matter for the statistics. To make the ids agree
everywhere, features with MS2 are placed first:

    feature_id 1..M      features with MS2 (exported to GNPS with the same ids)
    feature_id M+1..N    features without MS2 (internal tables only)
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pyopenms as oms

from atlas_ms.config import ExportSettings
from atlas_ms.preprocessing.msdata import (
    load_consensus_map,
    load_feature_map,
    set_parameters,
    store_consensus_map,
)

log = logging.getLogger(__name__)


def meta_value(element, key: str, default=None):
    """Meta value of an OpenMS object as a plain Python value (bytes decoded)."""
    if not element.metaValueExists(key):
        return default
    value = element.getMetaValue(key)
    return value.decode() if isinstance(value, bytes) else value


def order_features(consensus: oms.ConsensusMap, min_fraction: float) -> oms.ConsensusMap:
    """
    Drop features detected in too few runs, then put the features with MS2
    first (see module docstring). Within each group they are sorted by RT,
    then m/z.
    """
    n_runs = len(consensus.getColumnHeaders())
    kept = [cf for cf in consensus if cf.size() / n_runs >= min_fraction]
    kept.sort(key=lambda cf: (cf.getPeptideIdentifications().size() == 0, cf.getRT(), cf.getMZ()))
    ordered = oms.ConsensusMap(consensus)
    ordered.clear(False)  # empty the features, keep the column headers
    for cf in kept:
        ordered.push_back(cf)
    log.info("%d of %d features kept (min. detection fraction %.2f)", len(kept), consensus.size(), min_fraction)
    return ordered


def ms2_subset(consensus: oms.ConsensusMap) -> oms.ConsensusMap:
    """The features that have MS2 spectra (the first M features, see ``order_features``)."""
    subset = oms.ConsensusMap(consensus)
    subset.clear(False)
    for cf in consensus:
        if cf.getPeptideIdentifications().size() > 0:
            subset.push_back(cf)
    return subset


def gap_filled_features(feature_files: list[str | Path]) -> set[tuple[int, str]]:
    """
    (map index, feature id) of the features added by gap filling, read from
    the gap-filled feature maps (in map-index order). Consensus features only
    keep references to their features, not the features' meta values.
    """
    filled = set()
    for map_index, path in enumerate(feature_files):
        for feature in load_feature_map(path):
            if feature.metaValueExists("gap_filled"):
                filled.add((map_index, str(feature.getUniqueId())))
    return filled


def feature_tables(
    consensus: oms.ConsensusMap,
    run_names: list[str],
    gap_filled: set[tuple[int, str]] = frozenset(),
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build the features, quantification and gap-filled-flag tables from an
    ordered consensus map annotated with
    ``IonIdentityMolecularNetworking.annotateConsensusMap``.
    """
    n_features = consensus.size()
    intensities = np.full((n_features, len(run_names)), np.nan)
    filled = np.zeros((n_features, len(run_names)), dtype=bool)
    rows = []
    for row, cf in enumerate(consensus):
        for handle in cf.getFeatureList():
            intensities[row, handle.getMapIndex()] = handle.getIntensity()
            filled[row, handle.getMapIndex()] = (handle.getMapIndex(), str(handle.getUniqueId())) in gap_filled
        # Ion identity (IIMN): features that are adducts of the same molecule
        # list each other as "partners" (by row id) and share an "annotation
        # network number". The adduct itself ("best ion") is only evidence
        # when the feature has a partner. Alone, OpenMS just assumes [M+H]+
        # (or [M-H]-), so it is left empty.
        partners = str(meta_value(cf, "partners", "") or "")
        has_partner = bool(partners)
        rows.append({
            # annotateConsensusMap numbers rows 1..N in map order ("row ID").
            "feature_id": row + 1,
            "mz": cf.getMZ(),
            "rt": cf.getRT(),
            "charge": cf.getCharge(),
            "intensity": cf.getIntensity(),
            "quality": cf.getQuality(),
            "n_detected": cf.size(),  # runs with a value, detected or gap-filled
            "n_ms2": cf.getPeptideIdentifications().size(),
            "ion": meta_value(cf, "best ion") if has_partner else None,
            "ion_partners": partners.replace(",", ";") if has_partner else None,
            "ion_network": meta_value(cf, "annotation network number") if has_partner else None,
        })
    features = pd.DataFrame(rows)
    if features.empty:
        features = pd.DataFrame(columns=["feature_id", "mz", "rt", "charge", "intensity", "quality",
                                         "n_detected", "n_ms2", "ion", "ion_partners", "ion_network"])
    features["has_ms2"] = features["n_ms2"] > 0
    features["n_gap_filled"] = filled.sum(axis=1)
    quant = pd.DataFrame(intensities, columns=run_names)
    quant.insert(0, "feature_id", np.arange(1, n_features + 1))
    gap_filled_table = pd.DataFrame(filled, columns=run_names)
    gap_filled_table.insert(0, "feature_id", np.arange(1, n_features + 1))
    return features, quant, gap_filled_table


def gnps_metadata(samples: pd.DataFrame, mzml_files: list[str | Path]) -> pd.DataFrame:
    """
    Sample metadata in GNPS format: "filename" plus ATTRIBUTE_* columns.
    GNPS matches it to the quantification table by mzML file name.

    Written by its own rule (``gnps_metadata``): editing the metadata must
    not re-run the export, which would rewrite the MGF and re-run the
    spectral scoring downstream.
    """
    metadata = pd.DataFrame({"filename": [Path(f).name for f in mzml_files]})
    metadata["ATTRIBUTE_sample_type"] = samples["sample_type"].to_numpy()
    for column in samples.columns:
        if column.startswith("ATTRIBUTE_"):
            metadata[column] = samples[column].to_numpy()
    return metadata


def export_gnps(
    consensus: oms.ConsensusMap,
    mzml_files: list[str | Path],
    gnps_dir: str | Path,
    consensus_file: str | Path,
    settings: ExportSettings,
) -> None:
    """Write the GNPS FBMN files for the features with MS2 (except the metadata)."""
    gnps_dir = Path(gnps_dir)
    gnps_dir.mkdir(parents=True, exist_ok=True)
    subset = ms2_subset(consensus)
    # Row ids 1..M and IIMN adduct pairs, restricted to the exported features.
    oms.IonIdentityMolecularNetworking.annotateConsensusMap(subset)

    # GNPSMGFFile reads the consensus map from disk and takes each feature's
    # MS2 spectra from the mzML files (listed in map-index order).
    store_consensus_map(consensus_file, subset)
    mgf = oms.GNPSMGFFile()
    set_parameters(mgf, {"output_type": settings.ms2_spectrum})
    mgf.store(str(consensus_file), [str(f).encode() for f in mzml_files], str(gnps_dir / "ms2_spectra.mgf"))

    oms.GNPSQuantificationFile().store(subset, str(gnps_dir / "quantification_table.txt"))
    pairs = gnps_dir / "iimn_supplementary_pairs.csv"
    pairs.unlink(missing_ok=True)
    oms.IonIdentityMolecularNetworking.writeSupplementaryPairTable(subset, str(pairs))
    if not pairs.exists():
        # OpenMS writes nothing when no adduct pairs were found.
        pairs.write_text("ID1,ID2,EdgeType,Score,Annotation\n")

    # The MGF writer skips a feature when its chosen spectrum is empty (see
    # annotate.drop_empty_ms2, which prevents it). Report and check it.
    n_spectra = (gnps_dir / "ms2_spectra.mgf").read_text().count("BEGIN IONS")
    log.info("GNPS export: %d features with MS2, %d spectra written to the MGF", subset.size(), n_spectra)
    if n_spectra < subset.size():
        log.warning("%d features with MS2 have no spectrum in the MGF", subset.size() - n_spectra)


def export_results(
    consensus_file: str | Path,
    mzml_files: list[str | Path],
    run_names: list[str],
    features_out: str | Path,
    quant_out: str | Path,
    gap_filled_out: str | Path,
    gnps_dir: str | Path,
    gnps_consensus_out: str | Path,
    settings: ExportSettings,
    gap_filled_files: list[str | Path] | None = None,
) -> None:
    """
    File-level entry point. ``run_names`` are the sample names in map-index
    order; ``mzml_files`` and ``gap_filled_files`` (the gap-filled feature
    maps, if gap filling ran) list the runs' files in the same order.
    """
    consensus = order_features(load_consensus_map(consensus_file), settings.min_detection_fraction)

    # Internal tables: all features. Adduct partners are searched among all of
    # them, and row ids are the final feature ids.
    annotated = oms.ConsensusMap(consensus)
    oms.IonIdentityMolecularNetworking.annotateConsensusMap(annotated)
    gap_filled = gap_filled_features(gap_filled_files or [])
    features, quant, gap_filled_table = feature_tables(annotated, run_names, gap_filled)
    features.to_parquet(features_out, index=False)
    quant.to_parquet(quant_out, index=False)
    gap_filled_table.to_parquet(gap_filled_out, index=False)
    values = quant.drop(columns="feature_id").to_numpy()
    n_values = max(values.size, 1)
    log.info(
        "%d features (%d with MS2) x %d samples: %.1f%% missing values, %.1f%% gap-filled",
        len(features), int(features["has_ms2"].sum()), len(run_names),
        100 * np.isnan(values).sum() / n_values,
        100 * gap_filled_table.drop(columns="feature_id").to_numpy().sum() / n_values,
    )

    export_gnps(consensus, mzml_files, gnps_dir, gnps_consensus_out, settings)
