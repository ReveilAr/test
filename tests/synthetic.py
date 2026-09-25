"""
Synthetic LC-MS/MS data for tests.

Real raw files are large and can't be committed, so the tests generate
small centroided mzML files that look enough like a positive-mode lipidomics
DDA run to exercise the whole preprocessing chain:

* MS1 scans every 0.5 s, Gaussian elution profiles, M+1 / M+2 isotope peaks
  (so FeatureFinderMetabo sees real isotope patterns) and low-level noise;
* one MS2 scan at the apex of each compound (data-dependent acquisition),
  with class-specific fragments (e.g. m/z 184.0733 for phosphocholines);
* per-sample retention-time shifts (to exercise the alignment) and intensity
  changes (to give the statistics something to find);
* the option to make a compound too weak to be picked by the untargeted
  feature finder in one sample, so that gap filling has something to recover.

Everything is deterministic (fixed random seed) so tests are reproducible.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pyopenms as oms

# Mass difference between 13C and 12C: spacing of the isotope peaks.
C13_DELTA = 1.0033548


@dataclass
class Compound:
    """One analyte as seen in positive mode (m/z of the observed ion)."""

    name: str
    mz: float  # m/z of the monoisotopic ion
    rt: float  # apex retention time (s), before any per-sample shift
    intensity: float  # apex intensity of the monoisotopic peak
    fragments: list[tuple[float, float]] = field(default_factory=list)  # (m/z, relative intensity)
    isotopes: tuple[float, ...] = (1.0, 0.45, 0.12)  # relative M, M+1, M+2 (C40-C50 lipids)


# A handful of common lipids, with the ions and fragments expected in
# positive mode with an ammonium-containing mobile phase. PC 34:1 is present
# as both [M+H]+ and [M+Na]+ so that the adduct grouping can link them.
LIPIDS = [
    Compound("PC 34:1 [M+H]+", 760.5851, 180.0, 4e6, [(184.0733, 1.0), (104.1070, 0.05), (577.5190, 0.03)]),
    Compound("PC 34:1 [M+Na]+", 782.5670, 180.0, 8e5, [(723.4935, 1.0), (599.5010, 0.4), (146.9818, 0.2)]),
    Compound("SM 34:1;O2 [M+H]+", 703.5748, 150.0, 2e6, [(184.0733, 1.0), (264.2686, 0.08), (520.5088, 0.05)]),
    Compound("LPC 16:0 [M+H]+", 496.3398, 60.0, 3e6, [(184.0733, 1.0), (104.1070, 0.3), (478.3292, 0.1)]),
    Compound("Cer 34:1;O2 [M+H]+", 538.5194, 200.0, 1e6, [(264.2686, 1.0), (282.2791, 0.3), (520.5088, 0.4)]),
    Compound("TG 52:2 [M+NH4]+", 876.8014, 260.0, 5e6, [(577.5190, 1.0), (603.5347, 0.8), (859.7749, 0.05)]),
    Compound("CE 18:1 [M+NH4]+", 668.6344, 280.0, 1.5e6, [(369.3516, 1.0), (147.1168, 0.05)]),
    # Never selected for MS2 (no fragments): a feature without MS2.
    Compound("DG 34:1 [M+NH4]+", 612.5561, 230.0, 1.2e6, []),
]


def write_mzml(
    path: str | Path,
    compounds: list[Compound],
    rt_shift: float = 0.0,
    intensity_factors: dict[str, float] | None = None,
    seed: int = 0,
    rt_end: float = 320.0,
    scan_interval: float = 0.5,
    fwhm: float = 5.0,
    ppm_error: float = 1.0,
) -> Path:
    """
    Write one synthetic centroided DDA run and return its path.

    Parameters
    ----------
    compounds
        Analytes to simulate.
    rt_shift
        Seconds added to every apex: simulates a retention-time drift between runs.
    intensity_factors
        Optional multiplier per compound name (0 removes the compound from this run).
    seed
        Random seed (noise and mass errors), so each sample differs but reproducibly.
    fwhm
        Chromatographic peak width at half maximum (s).
    ppm_error
        Standard deviation of the random m/z error added to each centroid (ppm).
    """
    rng = np.random.default_rng(seed)
    factors = intensity_factors or {}
    sigma = fwhm / 2.3548  # Gaussian sigma from the FWHM

    # Apex position and height of every compound in this particular run.
    apexes = [(c, c.rt + rt_shift, c.intensity * factors.get(c.name, 1.0)) for c in compounds]
    apexes = [a for a in apexes if a[2] > 0]

    experiment = oms.MSExperiment()
    scan_rts = np.arange(0.0, rt_end, scan_interval)
    # DDA: each compound triggers one MS2 on the first MS1 scan at or after its apex.
    ms2_trigger = {c.name: scan_rts[np.searchsorted(scan_rts, apex_rt)] for c, apex_rt, _ in apexes}

    scan_number = 0
    for rt in scan_rts:
        # ---- MS1 scan: isotope envelopes of eluting compounds + noise ----
        mzs, intensities = [], []
        for compound, apex_rt, height in apexes:
            elution = height * np.exp(-0.5 * ((rt - apex_rt) / sigma) ** 2)
            if elution < 50:  # far from the apex: nothing measurable
                continue
            for k, rel in enumerate(compound.isotopes):
                mz = compound.mz + k * C13_DELTA
                mzs.append(mz * (1 + rng.normal(0, ppm_error) * 1e-6))
                intensities.append(elution * rel)
        # A few random low-intensity noise peaks, as seen in any real scan.
        noise_mz = rng.uniform(150, 1000, 15)
        mzs.extend(noise_mz)
        intensities.extend(rng.uniform(100, 600, noise_mz.size))

        scan_number += 1
        ms1 = oms.MSSpectrum()
        ms1.setMSLevel(1)
        ms1.setRT(float(rt))
        ms1.setNativeID(f"scan={scan_number}")
        ms1.setType(oms.SpectrumSettings.SpectrumType.CENTROID)
        order = np.argsort(mzs)  # peaks must be sorted by m/z
        ms1.set_peaks((np.asarray(mzs)[order], np.asarray(intensities, dtype=np.float32)[order]))
        experiment.addSpectrum(ms1)

        # ---- MS2 scans triggered on this MS1 scan ----
        n_ms2 = 0
        for compound, apex_rt, height in apexes:
            if ms2_trigger[compound.name] != rt or not compound.fragments:
                continue
            scan_number += 1
            n_ms2 += 1
            ms2 = oms.MSSpectrum()
            ms2.setMSLevel(2)
            # Consecutive MS2 scans are a few tens of ms apart, as on a real instrument.
            ms2.setRT(float(rt) + 0.05 * n_ms2)
            ms2.setNativeID(f"scan={scan_number}")
            ms2.setType(oms.SpectrumSettings.SpectrumType.CENTROID)
            precursor = oms.Precursor()
            # Instruments report the precursor with a small error, which the
            # pipeline's precursor correction step is meant to fix.
            precursor.setMZ(compound.mz * (1 + 3e-6))
            precursor.setCharge(1)
            precursor.setIntensity(height)
            ms2.setPrecursors([precursor])
            frag = sorted(compound.fragments)
            ms2.set_peaks((
                np.array([mz for mz, _ in frag]),
                np.array([rel * height * 0.2 for _, rel in frag], dtype=np.float32),
            ))
            experiment.addSpectrum(ms2)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    oms.MzMLFile().store(str(path), experiment)
    return path


def write_study(directory: str | Path) -> list[Path]:
    """
    Write a small 4-run study and return the mzML paths:

    * ``ctrl_1``, ``ctrl_2`` and ``treat_1``: the full lipid mix, with RT
      shifts of up to 4 s and a 3-fold TG increase in the treated sample;
    * ``treat_2``: same as ``treat_1`` but Cer is 300 times weaker, below
      the noise threshold used in the tests. Untargeted feature finding
      misses it there, so gap filling must recover it.
    """
    directory = Path(directory)
    runs = [
        ("ctrl_1", 0.0, {}, 1),
        ("ctrl_2", 2.0, {}, 2),
        ("treat_1", -3.0, {"TG 52:2 [M+NH4]+": 3.0}, 3),
        ("treat_2", 4.0, {"TG 52:2 [M+NH4]+": 3.0, "Cer 34:1;O2 [M+H]+": 1 / 300}, 4),
    ]
    return [
        write_mzml(directory / f"{name}.mzML", LIPIDS, rt_shift=shift, intensity_factors=factors, seed=seed)
        for name, shift, factors, seed in runs
    ]
