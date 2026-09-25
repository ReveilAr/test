"""Logging setup for the Snakemake rule scripts."""

import logging
import os
import subprocess
import sys
from pathlib import Path

import atlas_ms

# Keeps the log file open for the lifetime of the process (see log_to_file).
_log_file = None


def code_version() -> str:
    """
    Which code is running: package version, git commit (with "-dirty" if
    there are uncommitted changes) and the folder it is loaded from.

    Every rule log starts with this line, so a log always tells which version
    produced a result. The package is installed in editable mode, so it runs
    from whichever clone the conda environment was created from; that clone
    must be updated (`git pull`) to run new code.
    """
    folder = Path(atlas_ms.__file__).resolve().parent
    try:
        commit = subprocess.run(
            ["git", "describe", "--always", "--dirty"],
            cwd=folder, capture_output=True, text=True, timeout=10,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        commit = ""
    return f"ATLAS-MS {atlas_ms.__version__}, commit {commit or 'unknown'}, from {folder}"


def log_to_file(path: str) -> logging.Logger:
    """
    Send everything the current process prints to ``path``.

    OpenMS writes its progress messages from C++ straight to the process's
    standard output/error, bypassing Python. Redirecting the underlying file
    descriptors (1 and 2) captures those messages too, so each rule's log
    file tells the whole story and Snakemake's console stays readable.

    Snakemake runs every ``script:`` in its own Python process, so this never
    affects Snakemake itself.
    """
    global _log_file
    sys.stdout.flush()
    sys.stderr.flush()
    _log_file = open(path, "w", buffering=1)  # line-buffered: readable while the rule runs
    os.dup2(_log_file.fileno(), 1)
    os.dup2(_log_file.fileno(), 2)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr,
        force=True,
    )
    logger = logging.getLogger("atlas_ms")
    logger.info(code_version())
    return logger
