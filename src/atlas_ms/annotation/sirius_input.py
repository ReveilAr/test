"""
The SIRIUS input, prepared in the core environment.

SIRIUS runs in its own conda environment (``annotation/sirius.py``), which
has no pyOpenMS, so everything it needs is gathered here into one JSON file
(``work/sirius/input.json``), one entry per feature with an MS2 spectrum:

* ``feature_id``, precursor ``mz``, ``charge`` (signed), aligned ``rt_s``;
* ``adducts``: the ion found by the adduct grouping (e.g. "[M+NH4]+"),
  which SIRIUS then uses instead of trying every adduct;
* ``ms1``: the isotope pattern (m/z, intensity) from the run where the
  feature is most intense, at its apex. SIRIUS ranks molecular formulas
  with the isotope pattern and the fragmentation tree together, so the MS1
  pattern improves formula identification markedly;
* ``ms2``: the exported MS2 spectrum (the MGF one, as for the network).
"""

import json
import logging

import numpy as np
import pandas as pd

from atlas_ms.mgf import read_mgf
from atlas_ms.preprocessing.msdata import RunReader

log = logging.getLogger(__name__)


def sirius_input(features: pd.DataFrame, quant: pd.DataFrame, spectra: dict, reader: RunReader,
                 polarity: str, ppm: float) -> list[dict]:
    """One JSON-ready entry per feature with MS2 (see the module docstring)."""
    sign = 1 if polarity == "positive" else -1
    quant = quant.set_index("feature_id")
    entries = []
    for feature in features[features["feature_id"].isin(spectra)].itertuples():
        mz, intensity = spectra[feature.feature_id]
        charge = max(int(feature.charge), 1)
        intensities = quant.loc[feature.feature_id].dropna()
        pattern = []
        if not intensities.empty:
            sample = intensities.idxmax()  # the run where the feature is most intense
            pattern = reader.isotope_pattern(sample, feature.mz, feature.rt, ppm=ppm, charge=charge)
        ion = feature.ion if isinstance(feature.ion, str) else ""
        entries.append({
            "feature_id": int(feature.feature_id),
            "mz": float(feature.mz),
            "charge": sign * charge,
            "rt_s": float(feature.rt),
            "adducts": [ion] if ion else [],
            "ms1": [[float(a), float(b)] for a, b in pattern],
            "ms2": [[float(a), float(b)] for a, b in zip(mz, intensity)],
        })
    with_ms1 = sum(1 for entry in entries if len(entry["ms1"]) >= 2)
    log.info("SIRIUS input: %d features with MS2, %d with an isotope pattern (M and M+1 at least)",
             len(entries), with_ms1)
    return entries


def run_sirius_input(features_file, quant_file, mgf_file, samples: list[str], out, polarity: str, ppm: float) -> None:
    """File-level entry point of the ``prepare_sirius_input`` rule (run in the project folder)."""
    entries = sirius_input(pd.read_parquet(features_file), pd.read_parquet(quant_file), read_mgf(mgf_file),
                           RunReader(".", samples), polarity, ppm)
    with open(out, "w") as handle:
        json.dump(entries, handle, default=lambda value: value.item() if isinstance(value, np.generic) else value)
