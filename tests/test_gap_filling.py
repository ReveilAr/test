"""Gap-filling merge rules (pure logic, no raw data needed)."""

import pyopenms as oms
import pytest

from atlas_ms.preprocessing.gap_filling import (
    MIN_PAIRS_FOR_SCALE,
    intensity_scale,
    merge_gap_filled,
    monoisotopic_share,
)


def make_feature(mz: float, intensity: float, unique_id: int = 0, label: str | None = None) -> oms.Feature:
    feature = oms.Feature()
    feature.setMZ(mz)
    feature.setIntensity(intensity)
    if unique_id:
        feature.setUniqueId(unique_id)
    if label:
        feature.setMetaValue("label", label)
    return feature


def test_detected_values_are_kept_and_only_gaps_are_filled():
    features = oms.FeatureMap()
    features.push_back(make_feature(100.0, 1000.0, unique_id=1))  # complete consensus feature
    features.push_back(make_feature(200.0, 2000.0, unique_id=2))  # incomplete, detected in this run
    features.push_back(make_feature(999.0, 50.0, unique_id=9))    # in no consensus feature
    members = {"1": "", "2": "target_a"}

    extracted = oms.FeatureMap()
    extracted.push_back(make_feature(200.0, 1800.0, label="target_a"))  # detected here: not used
    extracted.push_back(make_feature(300.0, 500.0, label="target_b"))   # a gap: added

    merged, n_filled, scale = merge_gap_filled(features, members, extracted)
    values = sorted((round(f.getMZ()), f.getIntensity()) for f in merged)
    assert values == [(100, 1000.0), (200, 2000.0), (300, 500.0)]  # too few pairs: scale 1
    assert n_filled == 1 and scale == 1.0


def test_intensity_scale_is_the_median_ratio():
    # Detected values are 0.8 x the re-extracted ones, plus one outlier pair
    # (re-extraction took another peak): the median ignores it.
    n = MIN_PAIRS_FOR_SCALE + 5
    originals = {f"t{i}": 800.0 * (i + 1) for i in range(n)}
    extracted = {f"t{i}": 1000.0 * (i + 1) for i in range(n)}
    extracted["t0"] = 1.0
    scale, n_pairs = intensity_scale(originals, extracted)
    assert n_pairs == n
    assert scale == pytest.approx(0.8)


def test_monoisotopic_share():
    feature = oms.Feature()
    traces = []
    for isotope, probability in ((0, 0.7), (1, 0.3)):
        trace = oms.Feature()
        trace.setMetaValue("native_id", f"t_m500_z1_rt100_i{isotope}")
        trace.setMetaValue("isotope_probability", probability)
        traces.append(trace)
    feature.setSubordinates(traces)
    assert monoisotopic_share(feature) == pytest.approx(0.7)
    assert monoisotopic_share(oms.Feature()) == 1.0
