"""
Starting the Snakemake workflow on a project folder.

Used by the command line (`atlas-ms run`) and, later, by the app's Run button.
"""

import subprocess
import sys
from pathlib import Path

import atlas_ms

# The workflow sits next to the package in the repository. The package is
# installed in editable mode (see environment.yml), so it is found from here.
WORKFLOW_DIR = Path(atlas_ms.__file__).resolve().parents[2] / "workflow"

# Conda environments of the rules (e.g. ThermoRawFileParser) are created
# once per machine and shared by all projects.
CONDA_PREFIX = Path.home() / ".cache" / "atlas-ms" / "conda"


def snakemake_command(
    project_root: str | Path,
    cores: int = 4,
    dry_run: bool = False,
    extra_args: list[str] | None = None,
) -> list[str]:
    """The Snakemake command line that processes ``project_root``."""
    root = Path(project_root).resolve()
    command = [
        sys.executable, "-m", "snakemake",
        "--snakefile", str(WORKFLOW_DIR / "Snakefile"),
        "--directory", str(root),
        "--configfile", str(root / "project.yaml"),
        "--cores", str(cores),
        "--software-deployment-method", "conda",
        "--conda-prefix", str(CONDA_PREFIX),
        "--rerun-incomplete",  # redo steps interrupted by a crash or a Stop
    ]
    if dry_run:
        command.append("--dry-run")
    return command + list(extra_args or [])


def run_workflow(project_root: str | Path, cores: int = 4, dry_run: bool = False,
                 extra_args: list[str] | None = None) -> int:
    """Run the workflow in the foreground; returns Snakemake's exit code."""
    return subprocess.run(snakemake_command(project_root, cores, dry_run, extra_args)).returncode
