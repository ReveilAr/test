"""Merging rules of gap filling (pure logic, no raw data needed)."""

import pyopenms as oms

from atlas_ms.preprocessing.gap_filling import merge_gap_filled


def make_feature(mz: float, unique_id: int = 0, label: str | None = None) -> oms.Feature:
    feature = oms.Feature()
    feature.setMZ(mz)
    if unique_id:
        feature.setUniqueId(unique_id)
    if label:
        feature.setMetaValue("label", label)
    return feature


def test_original_feature_is_kept_when_re_extraction_fails():
    features = oms.FeatureMap()
    features.push_back(make_feature(100.0, unique_id=1))  # complete consensus feature
    features.push_back(make_feature(200.0, unique_id=2))  # incomplete, re-extracted below
    features.push_back(make_feature(300.0, unique_id=3))  # incomplete, re-extraction failed
    members = {"1": "", "2": "target_a", "3": "target_b"}

    extracted = oms.FeatureMap()
    extracted.push_back(make_feature(200.001, label="target_a"))  # found
    extracted.push_back(make_feature(400.0, label="target_c"))    # absent before, now filled

    merged = merge_gap_filled(features, members, extracted)
    # Kept: 100 (complete), 300 (safeguard); replaced: 200 -> 200.001; added: 400.
    assert sorted(round(f.getMZ(), 3) for f in merged) == [100.0, 200.001, 300.0, 400.0]
