"""
A small MGF reader, for code that must start fast (the app) or does not
need matchms (the lipid rules). matchms takes about 10 s to import.
"""

from pathlib import Path

import numpy as np


def read_mgf(path: str | Path) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """
    Read an MGF file into feature_id (the SCANS field) -> (m/z, intensity).
    The format is simple: BEGIN IONS, KEY=value lines, "m/z intensity" lines, END IONS.
    """
    spectra, scans, peaks = {}, None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line == "BEGIN IONS":
                scans, peaks = None, []
            elif line == "END IONS":
                if scans is not None:
                    array = np.array(peaks, dtype=float).reshape(-1, 2)
                    spectra[scans] = (array[:, 0], array[:, 1])
            elif line.startswith("SCANS="):
                scans = int(line.split("=", 1)[1])
            elif line and line[0].isdigit():
                mz, intensity = line.split()[:2]
                peaks.append((float(mz), float(intensity)))
    return spectra
