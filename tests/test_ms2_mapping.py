"""
MS2 spectra: attachment to features (instead of OpenMS's proteomics
IDMapper) and the spectrum written to the GNPS MGF.
"""

import numpy as np
import pyopenms as oms
import pytest

from atlas_ms.preprocessing.annotate import attach_ms2, ms2_spectra
from atlas_ms.preprocessing.export import binned_cosine, merge_spectra


def feature(mz: float, rt: float, charge: int = 1, rt_extent: tuple | None = None) -> oms.Feature:
    """A feature; with ``rt_extent``, a convex hull spanning those retention times."""
    f = oms.Feature()
    f.setMZ(mz)
    f.setRT(rt)
    f.setCharge(charge)
    if rt_extent:
        hull = oms.ConvexHull2D()
        hull.setHullPoints(np.array([[rt_extent[0], mz], [rt_extent[1], mz]], dtype=np.float32))
        f.setConvexHulls([hull])
    return f


def ms2(precursor_mz: float, rt: float, charge: int = 1, peaks=((184.0733, 100.0),)) -> oms.MSSpectrum:
    spectrum = oms.MSSpectrum()
    spectrum.setMSLevel(2)
    spectrum.setRT(rt)
    precursor = oms.Precursor()
    precursor.setMZ(precursor_mz)
    precursor.setCharge(charge)
    spectrum.setPrecursors([precursor])
    spectrum.set_peaks(([p[0] for p in peaks], [p[1] for p in peaks]))
    return spectrum


def test_each_spectrum_goes_to_the_feature_it_was_recorded_on():
    features = oms.FeatureMap()
    # Two co-eluting features 6.6 ppm apart (e.g. isobaric lipids), a third
    # one later, and one whose charge is unknown (0).
    features.push_back(feature(760.5851, 120.0, rt_extent=(110.0, 135.0)))
    features.push_back(feature(760.5901, 121.0, rt_extent=(112.0, 133.0)))
    features.push_back(feature(496.3398, 300.0))  # no hull: its apex +/- 5 s
    features.push_back(feature(703.5748, 200.0, charge=0))

    run = oms.MSExperiment()
    ms1 = oms.MSSpectrum()
    ms1.setMSLevel(1)
    run.addSpectrum(ms1)                          # 0: MS1, ignored
    run.addSpectrum(ms2(760.5895, 108.0))         # 1: closer to the second feature, 4 s before its hull
    run.addSpectrum(ms2(760.5851, 134.0))         # 2: first feature
    run.addSpectrum(ms2(760.5851, 150.0))         # 3: after both hulls (+5 s): unassigned
    run.addSpectrum(ms2(496.3400, 304.0, charge=0))  # 4: unknown precursor charge: kept
    run.addSpectrum(ms2(496.3398, 301.0, peaks=()))  # 5: no peaks: skipped
    run.addSpectrum(ms2(703.5750, 199.0))         # 6: unknown feature charge: kept
    run.addSpectrum(ms2(703.5750, 199.5, charge=2))  # 7: charge 2, but the feature's is unknown: kept
    attach_ms2(features, run)

    assert [ms2_spectra(f) for f in features] == [[2], [1], [4], [6, 7]]
    # Stored as a plain meta value, not as (proteomics) peptide identifications.
    assert all(f.getPeptideIdentifications().size() == 0 for f in features)


def test_charges_must_agree_when_both_are_known():
    features = oms.FeatureMap()
    features.push_back(feature(500.0, 100.0, charge=1))
    run = oms.MSExperiment()
    run.addSpectrum(ms2(500.0, 100.0, charge=2))
    attach_ms2(features, run)
    assert ms2_spectra(features[0]) == []


def test_merged_spectrum_keeps_the_shared_fragments():
    a = (np.array([184.0733, 264.2686, 520.5088]), np.array([100.0, 10.0, 5.0]))
    b = (np.array([184.0740, 264.2690, 400.0000]), np.array([1000.0, 80.0, 30.0]))  # 10x more intense run
    assert binned_cosine(a, a) == pytest.approx(1.0)
    assert binned_cosine(a, b) > 0.9
    mz, intensity = merge_spectra([a, b])
    # Each spectrum is scaled to its base peak, then averaged: fragments seen
    # once get half their relative intensity.
    assert mz == pytest.approx([184.07365, 264.2688, 400.0, 520.5088], abs=1e-3)
    assert intensity == pytest.approx([1.0, 0.09, 0.015, 0.025])
