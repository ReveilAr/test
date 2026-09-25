"""Project folders and sample tables."""

import pandas as pd
import pytest

from atlas_ms.project import Project, make_sample_table, validate_samples


def test_sample_names_are_safe_and_unique():
    table = make_sample_table(["/data/a/QC 01.mzML", "/data/b/QC 01.mzML", "/data/run.2.raw"])
    assert list(table["sample"]) == ["QC_01", "QC_01_2", "run_2"]
    assert set(table["sample_type"]) == {"sample"}


def test_unsupported_files_are_rejected():
    with pytest.raises(ValueError, match="Unsupported raw file"):
        make_sample_table(["/data/run.wiff"])


@pytest.mark.parametrize("change, message", [
    ({"sample_type": "qc"}, "Unknown sample_type"),
    ({"group": "a"}, "must start with 'ATTRIBUTE_'"),
    ({"sample": "bad name"}, "may only contain"),
])
def test_sample_table_validation(change, message):
    table = make_sample_table(["/data/run1.mzML"])
    for column, value in change.items():
        table[column] = value
    with pytest.raises(ValueError, match=message):
        validate_samples(table)


def test_create_project(tmp_path):
    raw = tmp_path / "run1.mzML"
    raw.touch()
    project = Project.create(tmp_path / "proj", [raw], instrument="qtof")
    assert project.load_config().presets.instrument == "qtof"
    samples = project.load_samples()
    assert list(samples.index) == ["run1"]
    assert samples.loc["run1", "file"] == str(raw)
    with pytest.raises(FileExistsError):
        Project.create(tmp_path / "proj", [raw])


def test_metadata_columns_survive_a_round_trip(tmp_path):
    raw = tmp_path / "run1.mzML"
    raw.touch()
    project = Project.create(tmp_path / "proj", [raw])
    samples = project.load_samples()
    samples["ATTRIBUTE_group"] = ["ctrl"]
    samples["sample_type"] = ["blank"]
    project.save_samples(samples)
    again = project.load_samples()
    pd.testing.assert_frame_equal(again, samples)
