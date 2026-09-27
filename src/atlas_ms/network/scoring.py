"""
Pairwise spectral similarity and candidate edges.

Two scores, chosen in the project parameters:

* modified cosine (matchms ``ModifiedCosineGreedy``): the GNPS score.
  Fragments are matched either directly or shifted by the precursor mass
  difference, so two molecules differing by a modification still match. It
  also gives the number of matched fragments. (Greedy fragment assignment,
  as in GNPS; matchms also offers an exact but slower Hungarian variant.)
* MS2DeepScore: a neural network trained on hundreds of thousands of
  library spectra predicts the structural similarity (Tanimoto) of the two
  molecules. Better at relating spectra that share few fragments; no
  matched-peak count.

Scoring every pair is the slow step, so it is done once and its result is
kept as a pool of *candidate* edges: for every spectrum, its
``candidates_per_spectrum`` best neighbours with at least
``min_candidate_score``. The network (``graph.py``) is then built from the
pool, and changing its cutoffs does not re-run the scoring.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from matchms import Spectrum

from atlas_ms.config import ScoringSettings
from atlas_ms.network.model_files import ms2deepscore_model_path

log = logging.getLogger(__name__)


def similarity_matrices(
    spectra: list[Spectrum], settings: ScoringSettings
) -> tuple[np.ndarray, np.ndarray | None]:
    """
    All-against-all scores as an N x N array, plus the N x N matched-fragment
    counts for the modified cosine (None for MS2DeepScore).
    """
    if settings.score == "modified_cosine":
        from matchms.similarity import ModifiedCosineGreedy

        result = ModifiedCosineGreedy(tolerance=settings.fragment_tolerance_da).matrix(
            spectra, spectra, is_symmetric=True, progress_bar=False
        )
        # A structured array with two fields: the score and the matched fragments.
        return result["score"].astype(np.float32), result["matches"].astype(np.int32)

    # Imported here: loading PyTorch takes a few seconds and is only needed now.
    from ms2deepscore import MS2DeepScore
    from ms2deepscore.models import load_model

    # allow_legacy: the published pretrained model uses an older file format.
    model = load_model(ms2deepscore_model_path(settings), allow_legacy=True)
    scores = MS2DeepScore(model, progress_bar=False).matrix(spectra, spectra, is_symmetric=True, progress_bar=False)
    return np.asarray(scores, dtype=np.float32), None


def candidate_edges(
    scores: np.ndarray,
    matches: np.ndarray | None,
    feature_ids: list[int],
    settings: ScoringSettings,
) -> pd.DataFrame:
    """
    The pool of candidate edges from a full score matrix: for each spectrum,
    its best ``candidates_per_spectrum`` neighbours scoring at least
    ``min_candidate_score``. An edge found from both ends is listed once.

    Columns: source, target (feature ids, source < target), score, matched_peaks.
    """
    n = scores.shape[0]
    if n < 2:
        return pd.DataFrame(columns=["source", "target", "score", "matched_peaks"])
    scores = scores.copy()
    np.fill_diagonal(scores, -1.0)  # a spectrum is not its own neighbour
    k = min(settings.candidates_per_spectrum, n - 1)
    # Indices of the k best neighbours of every row, without a full sort.
    best = np.argpartition(-scores, k - 1, axis=1)[:, :k]
    rows = np.repeat(np.arange(n), k)
    cols = best.ravel()
    keep = scores[rows, cols] >= settings.min_candidate_score
    rows, cols = rows[keep], cols[keep]
    # Undirected edges: order each pair, then drop duplicates.
    first, second = np.minimum(rows, cols), np.maximum(rows, cols)
    pairs = np.unique(np.stack([first, second], axis=1), axis=0)
    ids = np.asarray(feature_ids)
    return pd.DataFrame({
        "source": ids[pairs[:, 0]],
        "target": ids[pairs[:, 1]],
        "score": scores[pairs[:, 0], pairs[:, 1]],
        "matched_peaks": matches[pairs[:, 0], pairs[:, 1]] if matches is not None else np.nan,
    })


def score_spectra(spectra: list[Spectrum], settings: ScoringSettings) -> pd.DataFrame:
    """Candidate edges between the QC-passed spectra of a project."""
    feature_ids = [spectrum.get("feature_id") for spectrum in spectra]
    if len(spectra) < 2:
        return candidate_edges(np.zeros((len(spectra), len(spectra))), None, feature_ids, settings)
    scores, matches = similarity_matrices(spectra, settings)
    edges = candidate_edges(scores, matches, feature_ids, settings)
    log.info("%s: %d spectra, %d candidate edges", settings.score, len(spectra), len(edges))
    return edges


def run_scoring(spectra_file: str | Path, candidates_out: str | Path, settings: ScoringSettings,
                threads: int = 1) -> None:
    """File-level entry point: candidate edges of the QC-passed spectra."""
    from atlas_ms.network.spectra import load_spectra

    if settings.score == "ms2deepscore":
        import torch

        torch.set_num_threads(threads)  # the cores Snakemake gave this step

    score_spectra(load_spectra(spectra_file), settings).to_parquet(candidates_out, index=False)
