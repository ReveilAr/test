"""
Retention-time alignment of all runs (UmetaFlow's `MapAligner` rule).

Retention times drift a little between injections. MapAlignerPoseClustering
pairs the features of each run with those of a reference run (similar m/z
and RT) and fits a linear transformation of the RT axis. The run with the
most features is the reference.

Only the transformations are saved (one small trafoXML per run). The steps
that need aligned data apply them on the fly (see ``msdata.load_run``).
"""

import logging
from pathlib import Path

import numpy as np
import pyopenms as oms

from atlas_ms.config import AlignmentSettings
from atlas_ms.preprocessing.msdata import load_feature_map, set_parameters, store_trafo

log = logging.getLogger(__name__)


def identity_trafo() -> oms.TransformationDescription:
    """A transformation that leaves retention times unchanged."""
    trafo = oms.TransformationDescription()
    trafo.fitModel("identity")
    return trafo


def align_runs(
    feature_files: list[str | Path],
    trafo_files: list[str | Path],
    settings: AlignmentSettings,
) -> None:
    """Compute one RT transformation per run (same order as ``feature_files``)."""
    feature_maps = [load_feature_map(path) for path in feature_files]
    sizes = [fmap.size() for fmap in feature_maps]
    reference = int(np.argmax(sizes))
    log.info("Reference run: %s (%d features)", feature_files[reference], sizes[reference])

    aligner = oms.MapAlignmentAlgorithmPoseClustering()
    set_parameters(aligner, {
        "max_num_peaks_considered": -1,  # use all features, not only the most intense 1000
        "superimposer:mz_pair_max_distance": 0.05,  # Da, as in UmetaFlow
        "pairfinder:distance_MZ:max_difference": settings.mz_max_ppm,
        "pairfinder:distance_MZ:unit": "ppm",
        "pairfinder:distance_RT:max_difference": settings.rt_max_difference_s,
    })
    aligner.setReference(feature_maps[reference])

    for index, (fmap, trafo_file) in enumerate(zip(feature_maps, trafo_files)):
        if index == reference or fmap.size() == 0:
            # The reference defines the time axis. An empty run (e.g. a clean
            # blank) has nothing to align and OpenMS would refuse it.
            if fmap.size() == 0:
                log.warning("%s has no features: not aligned", feature_files[index])
            trafo = identity_trafo()
        else:
            trafo = oms.TransformationDescription()
            aligner.align(fmap, trafo)
            log.info(
                "%s: %s RT model from %d feature pairs",
                feature_files[index], trafo.getModelType(), len(trafo.getDataPoints()),
            )
        store_trafo(trafo_file, trafo)
