"""
Starting the Snakemake workflow on a project folder.

Used by the command line (`atlas-ms run`) and by the app's Run button.

Parallel work. ``--cores N`` is a *maximum*: Snakemake runs as many jobs at
once as fit in N cores (one core per job, unless a rule asks for more), so a
step with a single job uses one core, and 20 runs to process use up to N.
Two things would make more parallel jobs slower rather than faster, and are
prevented here:

* memory: steps that load a whole run (feature finding, MS2 mapping, gap
  filling) declare an estimate of their memory use (``mem_mb``), and
  Snakemake is given 80% of the memory currently available, so it starts
  fewer runs at once on a laptop instead of swapping;
* thread oversubscription: libraries such as OpenMP or BLAS start one thread
  per core in *every* job. They are limited to one thread per job, except
  where a rule asks for more: the steps that process one run share the
  cores between the runs (4 runs on 16 cores: 4 OpenMP threads each), and
  the steps that are parallel inside (MS2DeepScore, SIRIUS, MS2Query) take
  all the cores for themselves (``threads`` in their rules).
"""

import os
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


def available_cores() -> int:
    """CPU cores this process may use (respects limits set with taskset or by a cluster)."""
    return len(os.sched_getaffinity(0))


def available_memory_mb() -> int | None:
    """
    Memory (MB) that can be used without swapping: MemAvailable of
    /proc/meminfo (Linux). None if unknown.
    """
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) // 1024  # the file gives kB
    except OSError:
        pass
    return None


def snakemake_env() -> dict[str, str]:
    """
    The environment of the Snakemake process (and of its jobs): one OpenMP /
    BLAS thread per job (see the module docstring). Rules that are parallel
    inside set their own thread count.
    """
    env = dict(os.environ)
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[variable] = "1"
    return env


def snakemake_command(
    project_root: str | Path,
    cores: int | None = None,
    dry_run: bool = False,
    extra_args: list[str] | None = None,
) -> list[str]:
    """
    The Snakemake command line that processes ``project_root``, with at
    most ``cores`` cores (default: all available) and 80% of the available
    memory.
    """
    root = Path(project_root).resolve()
    command = [
        sys.executable, "-m", "snakemake",
        "--snakefile", str(WORKFLOW_DIR / "Snakefile"),
        "--directory", str(root),
        "--configfile", str(root / "project.yaml"),
        "--cores", str(cores or available_cores()),
        "--software-deployment-method", "conda",
        "--conda-prefix", str(CONDA_PREFIX),
        "--rerun-incomplete",  # redo steps interrupted by a crash or a Stop
    ]
    memory = available_memory_mb()
    if memory:
        command += ["--resources", f"mem_mb={int(0.8 * memory)}"]
    if dry_run:
        command.append("--dry-run")
    return command + list(extra_args or [])


def run_workflow(project_root: str | Path, cores: int | None = None, dry_run: bool = False,
                 extra_args: list[str] | None = None) -> int:
    """Run the workflow in the foreground; returns Snakemake's exit code."""
    command = snakemake_command(project_root, cores, dry_run, extra_args)
    return subprocess.run(command, env=snakemake_env()).returncode
