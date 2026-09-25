"""Shared test fixtures: a synthetic study, processed once per test session."""

import subprocess
import sys
from pathlib import Path

import pytest

from atlas_ms.project import Project
from atlas_ms.runner import WORKFLOW_DIR

sys.path.insert(0, str(Path(__file__).parent))
from synthetic import write_study  # noqa: E402


def snakemake(project_root: Path, *args: str) -> subprocess.CompletedProcess:
    """
    Run the workflow on a project, without conda. The synthetic runs are
    mzML files, so the only conda rule (Thermo conversion) is never needed.
    """
    command = [
        sys.executable, "-m", "snakemake",
        "--snakefile", str(WORKFLOW_DIR / "Snakefile"),
        "--directory", str(project_root),
        "--configfile", str(project_root / "project.yaml"),
        "--cores", "2",
        *args,
    ]
    return subprocess.run(command, capture_output=True, text=True)


@pytest.fixture(scope="session")
def study_files(tmp_path_factory) -> list[Path]:
    """Four synthetic mzML runs (see synthetic.write_study)."""
    return write_study(tmp_path_factory.mktemp("raw"))


@pytest.fixture(scope="session")
def processed_project(tmp_path_factory, study_files) -> Project:
    """A project made from the synthetic study, fully processed by the workflow."""
    project = Project.create(tmp_path_factory.mktemp("projects") / "study", study_files)
    # The synthetic spectra have 2-3 fragments: relax the fragment-count
    # cutoffs, which are meant for real spectra. With a single matched
    # fragment, TG and CE would be linked by coincidence: CE's m/z 369.35
    # shifted by their precursor difference (208.17) lands on TG's m/z 577.52,
    # giving a modified cosine of 0.78. That is why GNPS asks for several
    # matched fragments.
    config = project.load_config()
    config.spectrum_qc.min_peaks = 2
    config.network.min_matched_peaks = 2
    project.save_config(config)
    samples = project.load_samples()
    samples["ATTRIBUTE_group"] = ["ctrl", "ctrl", "treat", "treat"]
    project.save_samples(samples)
    result = snakemake(project.root)
    assert result.returncode == 0, result.stderr[-5000:]
    return project
