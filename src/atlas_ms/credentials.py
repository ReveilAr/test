"""
The SIRIUS account, kept outside every project.

The app's Setup tab (SIRIUS card) writes it; the ``run_sirius`` rule reads it
and logs SIRIUS in with it when SIRIUS is not logged in yet. It lives in
``~/.config/atlas-ms/sirius_account.json``, readable by you only (like the
credential files of many command-line tools), and never in project.yaml or
in the rule parameters (Snakemake keeps those in the project folder): a
project folder can be shared without the password.

After the first login, SIRIUS keeps its own session, so the password is only
used again if that session ends.
"""

import json
import os
from pathlib import Path

ACCOUNT_FILE = Path.home() / ".config" / "atlas-ms" / "sirius_account.json"


def save_sirius_account(username: str, password: str, accept_terms: bool, path: Path = ACCOUNT_FILE) -> None:
    """Store the account; the file is created readable and writable by its owner only."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.unlink(missing_ok=True)  # recreated below with owner-only permissions from the start
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(descriptor, "w") as handle:
        json.dump({"username": username, "password": password, "accept_terms": bool(accept_terms)}, handle)


def load_sirius_account(path: Path = ACCOUNT_FILE) -> dict | None:
    """The stored account, or None."""
    path = Path(path)
    return json.loads(path.read_text()) if path.is_file() else None


def forget_sirius_account(path: Path = ACCOUNT_FILE) -> None:
    Path(path).unlink(missing_ok=True)
