"""
Where the MS2DeepScore model file lives, and its download.

Kept apart from ``scoring.py`` on purpose: the Snakefile needs these paths
when it starts, and importing matchms (as ``scoring.py`` does) takes about
10 s, which every Snakemake start (even a dry run) would pay.
"""

import logging
import urllib.request
from pathlib import Path

from atlas_ms.config import ScoringSettings

log = logging.getLogger(__name__)

# Pretrained MS2DeepScore model (both ion modes; see the MS2DeepScore README).
MS2DEEPSCORE_MODEL_URL = "https://zenodo.org/records/17826815/files/ms2deepscore_model.pt?download=1"
DEFAULT_MODEL_PATH = Path.home() / ".cache" / "atlas-ms" / "models" / "ms2deepscore_model.pt"


def ms2deepscore_model_path(settings: ScoringSettings) -> Path:
    """The model file to use: the one given in the parameters, or the pretrained one."""
    return Path(settings.ms2deepscore_model).expanduser() if settings.ms2deepscore_model else DEFAULT_MODEL_PATH


def download_model(path: str | Path, url: str = MS2DEEPSCORE_MODEL_URL) -> None:
    """Download the pretrained MS2DeepScore model (about 130 MB), once."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(".part")  # a failed download never looks complete
    log.info("Downloading %s", url)
    urllib.request.urlretrieve(url, partial)
    partial.rename(path)
