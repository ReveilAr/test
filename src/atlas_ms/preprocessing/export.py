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

MS2 spectra: each feature's spectra were attached in each run
(``annotate.attach_ms2``). The consensus map only refers to its features
(map index + feature id), so the linked feature maps are read again to
collect them (``linked_features``). ``write_mgf`` then writes one spectrum
per feature in the MGF format of GNPS FBMN, with the fields of OpenMS's
GNPSMGFFile (which UmetaFlow uses, and which needs the MS2 spectra stored as
proteomics "peptide identifications").

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
from atlas_ms.preprocessing.annotate import ms2_spectra
from atlas_ms.preprocessing.msdata import load_consensus_map, load_feature_map

log = logging.getLogger(__name__)


def meta_value(element, key: str, default=None):
    """Meta value of an OpenMS object as a plain Python value (bytes decoded)."""
    if not element.metaValueExists(key):
        return default
    value = element.getMetaValue(key)
    return value.decode() if isinstance(value, bytes) else value


def linked_features(feature_files: list[str | Path]) -> tuple[dict, set]:
    """
    What the consensus map does not keep about its features, read from the
    linked feature maps (in map-index order):

    * the MS2 spectra attached to each feature:
      (map index, feature id) -> spectrum indices in that run;
    * the features added by gap filling: a set of (map index, feature id).
    """
    ms2, gap_filled = {}, set()
    for map_index, path in enumerate(feature_files):
        for feature in load_feature_map(path):
            key = (map_index, str(feature.getUniqueId()))
            spectra = ms2_spectra(feature)
            if spectra:
                ms2[key] = spectra
            if feature.metaValueExists("gap_filled"):
                gap_filled.add(key)
    return ms2, gap_filled


def consensus_ms2(consensus: oms.ConsensusMap, ms2: dict) -> list[list[tuple[int, int]]]:
    """
    The MS2 spectra of every consensus feature, as (map index, spectrum
    index) pairs: those of its features in every run (``linked_features``).
    """
    return [[(handle.getMapIndex(), index)
             for handle in cf.getFeatureList()
             for index in ms2.get((handle.getMapIndex(), str(handle.getUniqueId())), [])]
            for cf in consensus]


def order_features(consensus: oms.ConsensusMap, ms2: list, min_fraction: float) -> tuple[oms.ConsensusMap, list]:
    """
    Drop features detected in too few runs, then put the features with MS2
    first (see module docstring). Within each group they are sorted by RT,
    then m/z. ``ms2`` (one list per consensus feature, ``consensus_ms2``)
    is returned in the new order.
    """
    n_runs = len(consensus.getColumnHeaders())
    features = list(consensus)  # copies (pyOpenMS), made once
    kept = [i for i, cf in enumerate(features) if cf.size() / n_runs >= min_fraction]
    kept.sort(key=lambda i: (not ms2[i], features[i].getRT(), features[i].getMZ()))
    ordered = oms.ConsensusMap(consensus)
    ordered.clear(False)  # empty the features, keep the column headers
    for i in kept:
        ordered.push_back(features[i])
    log.info("%d of %d features kept (min. detection fraction %.2f)", len(kept), consensus.size(), min_fraction)
    return ordered, [ms2[i] for i in kept]


def ms2_subset(consensus: oms.ConsensusMap, ms2: list) -> oms.ConsensusMap:
    """The features that have MS2 spectra (the first M features, see ``order_features``)."""
    subset = oms.ConsensusMap(consensus)
    subset.clear(False)
    for cf, spectra in zip(consensus, ms2):
        if spectra:
            subset.push_back(cf)
    return subset


def feature_tables(
    consensus: oms.ConsensusMap,
    run_names: list[str],
    ms2: list,
    gap_filled: set[tuple[int, str]] = frozenset(),
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build the features, quantification and gap-filled-flag tables from an
    ordered consensus map annotated with
    ``IonIdentityMolecularNetworking.annotateConsensusMap``, and its MS2
    spectra (one list per feature, ``order_features``).
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
            "n_ms2": len(ms2[row]),
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


# ---- MGF of the MS2 spectra ----------------------------------------------------

# Merging of similar spectra ("merged_spectra"), with OpenMS GNPSMGFFile's defaults.
MERGE_MAX_RUNS = 5  # best spectrum of at most this many runs (the most intense ones)
MERGE_MIN_COSINE = 0.9  # similarity to the best spectrum needed to be merged
MERGE_BIN_DA = 0.02  # fragments closer than this are the same fragment


class RunSpectra:
    """
    MS2 spectra read from the runs' mzML files on demand. The files are
    opened "on disc" (OpenMS reads one spectrum through the file index), so a
    100-run study is never loaded in memory.
    """

    def __init__(self, mzml_files: list[str | Path]):
        self.mzml_files = [str(f) for f in mzml_files]
        self.runs = {}

    def get(self, map_index: int, spectrum_index: int) -> tuple[np.ndarray, np.ndarray, float]:
        """(m/z, intensities, retention time in s) of one spectrum."""
        if map_index not in self.runs:
            run = oms.OnDiscMSExperiment()
            if not run.openFile(self.mzml_files[map_index]):
                raise OSError(f"Cannot read the spectra of {self.mzml_files[map_index]} (indexed mzML needed)")
            self.runs[map_index] = run
        spectrum = self.runs[map_index].getSpectrum(spectrum_index)
        mz, intensity = spectrum.get_peaks()
        return mz, intensity, spectrum.getRT()


def _bins(mz: np.ndarray, intensity: np.ndarray) -> dict[int, float]:
    """A spectrum's intensity summed in bins of ``MERGE_BIN_DA``."""
    keys, inverse = np.unique(np.floor(mz / MERGE_BIN_DA).astype(np.int64), return_inverse=True)
    return dict(zip(keys.tolist(), np.bincount(inverse, weights=intensity).tolist()))


def binned_cosine(a: tuple, b: tuple) -> float:
    """Cosine similarity of two spectra (m/z, intensities) on ``MERGE_BIN_DA`` bins."""
    bins_a, bins_b = _bins(*a), _bins(*b)
    dot = sum(value * bins_b.get(key, 0.0) for key, value in bins_a.items())
    norm = np.sqrt(sum(v * v for v in bins_a.values()) * sum(v * v for v in bins_b.values()))
    return dot / norm if norm > 0 else 0.0


def merge_spectra(spectra: list[tuple]) -> tuple[np.ndarray, np.ndarray]:
    """
    One spectrum from several spectra of the same compound: each is scaled
    to its base peak, then fragments closer than ``MERGE_BIN_DA`` are merged
    (intensity-weighted m/z, intensity averaged over the spectra: a fragment
    seen in one spectrum only is weakened, which removes noise).
    """
    mz = np.concatenate([s[0] for s in spectra])
    intensity = np.concatenate([s[1] / s[1].max() for s in spectra])
    order = np.argsort(mz)
    mz, intensity = mz[order], intensity[order]
    # A new fragment starts wherever the gap to the previous peak exceeds the bin width.
    group = np.cumsum(np.diff(mz, prepend=-np.inf) > MERGE_BIN_DA) - 1
    total = np.bincount(group, weights=intensity)
    return np.bincount(group, weights=mz * intensity) / total, total / len(spectra)


def feature_spectrum(cf: oms.ConsensusFeature, refs: list[tuple[int, int]], runs: RunSpectra,
                     mode: str) -> tuple[np.ndarray, np.ndarray, int, float]:
    """
    The MS2 spectrum exported for one feature: (m/z, intensities, spectrum
    index, retention time) of the best spectrum, merged with others in
    "merged_spectra" mode.

    The best spectrum comes from the run where the feature is most intense
    (among the runs that have MS2 spectra for it): the strongest precursor
    gives the cleanest fragments. Within that run, it is the spectrum with
    the highest total intensity (runs usually have one or two per feature).
    """
    intensity_in = {handle.getMapIndex(): handle.getIntensity() for handle in cf.getFeatureList()}
    by_run = {}  # map index -> spectra of the feature in that run
    for map_index, index in refs:
        by_run.setdefault(map_index, []).append((index, *runs.get(map_index, index)))
    ordered_runs = sorted(by_run, key=lambda m: -intensity_in.get(m, 0.0))
    # Best spectrum of each run: highest total fragment intensity.
    best = [max(by_run[m], key=lambda s: s[2].sum()) for m in ordered_runs]
    index, mz, intensity, rt = best[0]
    if mode == "merged_spectra":
        similar = [(s[1], s[2]) for s in best[:MERGE_MAX_RUNS]
                   if binned_cosine((mz, intensity), (s[1], s[2])) >= MERGE_MIN_COSINE]
        mz, intensity = merge_spectra(similar)  # the best spectrum is always in (cosine 1)
    return mz, intensity, index, rt


def write_mgf(consensus: oms.ConsensusMap, ms2: list, mzml_files: list[str | Path], path: str | Path,
              mode: str) -> int:
    """
    The MS2 spectra of the exported features (row ids 1..M in map order),
    one per feature, in the MGF format of GNPS FBMN. The fields are those of
    OpenMS's GNPSMGFFile ("OpenMS" input of GNPS): SCANS = row id, PEPMASS =
    the feature's m/z, FILE_INDEX and RTINSECONDS = the spectrum's index and
    retention time in its mzML file. Returns the number of spectra written.
    """
    runs = RunSpectra(mzml_files)
    n_written = 0
    with open(path, "w") as out:
        for row, (cf, refs) in enumerate(zip(consensus, ms2), start=1):
            mz, intensity, index, rt = feature_spectrum(cf, refs, runs, mode)
            keep = intensity > 0
            # Charge: the highest of the runs (unknown = 1), sign as in OpenMS.
            charge = max([cf.getCharge(), *(h.getCharge() for h in cf.getFeatureList())])
            lines = [
                "BEGIN IONS",
                f"OUTPUT={mode}",
                f"SCANS={row}",
                f"FEATURE_ID=e_{cf.getUniqueId()}",
                "MSLEVEL=2",
                f"CHARGE={abs(charge) or 1}{'+' if charge >= 0 else '-'}",
                f"PEPMASS={cf.getMZ():.6f}",
                f"FILE_INDEX={index}",
                f"RTINSECONDS={rt:.3f}",
                *(f"{m:.4f}\t{i:.4f}" for m, i in zip(mz[keep], intensity[keep])),
                "END IONS",
                "",
            ]
            out.write("\n".join(lines) + "\n")
            n_written += 1
    return n_written


def export_gnps(
    consensus: oms.ConsensusMap,
    ms2: list,
    mzml_files: list[str | Path],
    gnps_dir: str | Path,
    settings: ExportSettings,
) -> None:
    """
    Write the GNPS FBMN files for the features with MS2 (except the
    metadata). ``ms2`` lists each feature's spectra (``order_features``).
    """
    gnps_dir = Path(gnps_dir)
    gnps_dir.mkdir(parents=True, exist_ok=True)
    subset = ms2_subset(consensus, ms2)
    # Row ids 1..M and IIMN adduct pairs, restricted to the exported features.
    oms.IonIdentityMolecularNetworking.annotateConsensusMap(subset)

    n_spectra = write_mgf(subset, [refs for refs in ms2 if refs], mzml_files, gnps_dir / "ms2_spectra.mgf",
                          settings.ms2_spectrum)
    oms.GNPSQuantificationFile().store(subset, str(gnps_dir / "quantification_table.txt"))
    pairs = gnps_dir / "iimn_supplementary_pairs.csv"
    pairs.unlink(missing_ok=True)
    oms.IonIdentityMolecularNetworking.writeSupplementaryPairTable(subset, str(pairs))
    if not pairs.exists():
        # OpenMS writes nothing when no adduct pairs were found.
        pairs.write_text("ID1,ID2,EdgeType,Score,Annotation\n")
    log.info("GNPS export: %d features with MS2, %d spectra written to the MGF", subset.size(), n_spectra)


def export_results(
    consensus_file: str | Path,
    feature_files: list[str | Path],
    mzml_files: list[str | Path],
    run_names: list[str],
    features_out: str | Path,
    quant_out: str | Path,
    gap_filled_out: str | Path,
    gnps_dir: str | Path,
    settings: ExportSettings,
) -> None:
    """
    File-level entry point. ``run_names`` are the sample names in map-index
    order; ``feature_files`` (the feature maps that were linked into
    ``consensus_file``) and ``mzml_files`` list the runs' files in the same
    order.
    """
    consensus = load_consensus_map(consensus_file)
    ms2, gap_filled = linked_features(feature_files)
    if not ms2 and any(cf.getPeptideIdentifications().size() for cf in consensus):
        # Feature maps from before attach_ms2: the MS2 spectra would be lost.
        raise RuntimeError("These feature maps come from an older ATLAS-MS version (MS2 spectra stored as "
                           "OpenMS peptide identifications). Redo the MS2 mapping: "
                           "atlas-ms run <project> -- --forcerun annotate_run")
    consensus, ms2 = order_features(consensus, consensus_ms2(consensus, ms2), settings.min_detection_fraction)

    # Internal tables: all features. Adduct partners are searched among all of
    # them, and row ids are the final feature ids.
    annotated = oms.ConsensusMap(consensus)
    oms.IonIdentityMolecularNetworking.annotateConsensusMap(annotated)
    features, quant, gap_filled_table = feature_tables(annotated, run_names, ms2, gap_filled)
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

    export_gnps(consensus, ms2, mzml_files, gnps_dir, settings)
