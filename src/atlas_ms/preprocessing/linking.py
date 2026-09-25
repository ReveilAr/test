"""
Linking features across runs (UmetaFlow's `FeatureLinker` rules).

FeatureLinkerUnlabeledKD groups the features of all runs that have the same
m/z and (aligned) RT into "consensus features": one row of the final feature
table, with one intensity per sample. The KD-tree variant is fast even for
large studies.

The order of the runs defines the "map index" used by OpenMS in every
multi-run file. It is always the order of ``samples.tsv``.
"""

import logging
from pathlib import Path

import pyopenms as oms

from atlas_ms.config import LinkingSettings
from atlas_ms.preprocessing.msdata import load_feature_map, set_parameters, store_consensus_map

log = logging.getLogger(__name__)


def link_feature_maps(
    feature_maps: list[oms.FeatureMap],
    run_names: list[str],
    mzml_files: list[str | Path],
    settings: LinkingSettings,
) -> oms.ConsensusMap:
    """Group corresponding features of all runs into one ConsensusMap."""
    consensus = oms.ConsensusMap()

    # One column per run. The mzML path is what the GNPS export reads the MS2
    # spectra from, and its file name is how GNPS matches the sample metadata.
    headers = consensus.getColumnHeaders()
    for index, (fmap, name, mzml) in enumerate(zip(feature_maps, run_names, mzml_files)):
        header = headers.get(index, oms.ColumnHeader())
        header.filename = str(mzml)
        header.label = name
        header.size = fmap.size()
        header.unique_id = fmap.getUniqueId()
        headers[index] = header
    consensus.setColumnHeaders(headers)

    linker = oms.FeatureGroupingAlgorithmKD()
    set_parameters(linker, {
        "warp:enabled": "false",  # the runs are already aligned
        "link:rt_tol": settings.rt_tol_s,
        "link:mz_tol": settings.mz_tol_ppm,
        "mz_unit": "ppm",
    })
    linker.group(feature_maps, consensus)
    consensus.setUniqueIds()
    log.info("%d runs linked into %d consensus features", len(feature_maps), consensus.size())
    return consensus


def link_runs(
    feature_files: list[str | Path],
    run_names: list[str],
    mzml_files: list[str | Path],
    consensus_out: str | Path,
    settings: LinkingSettings,
) -> None:
    """File-level entry point of ``link_feature_maps``."""
    feature_maps = [load_feature_map(path) for path in feature_files]
    store_consensus_map(consensus_out, link_feature_maps(feature_maps, run_names, mzml_files, settings))
