"""Logging setup for the Snakemake rule scripts."""

import logging
import os
import sys

# Keeps the log file open for the lifetime of the process (see log_to_file).
_log_file = None


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
    return logging.getLogger("atlas_ms")
