"""
Reading spectral library files, one spectrum at a time.

Formats:

* MSP (NIST, MS-DIAL, MoNA...) and MGF (GNPS...): read by matchms.
* JSON: three different formats share the extension, so each record is
  recognised by its keys:

  - GNPS (``peaks_json``): converted by matchms;
  - MoNA exports (``compound`` + ``spectrum``, e.g. ``MoNA-export-LipidBlast.json``);
  - MassBank (``PK$PEAK``, the ``MassBank.json`` of the MassBank-data
    releases, keys of the MassBank record format).

  matchms only reads the GNPS format, and a MoNA or MassBank file gave no
  spectrum at all (without an error). They are converted here.

JSON files are streamed (``json_records``): ``json.load`` builds every record
in memory at once, several times the size of the file, and MoNA exports are
gigabytes.

Retention times of MoNA and MassBank records are converted to minutes
(their records say "5.68 min" or "340 s"), so keep ``rt_unit`` = "min" for
these libraries.
"""

import json
import logging
import re
from collections.abc import Iterator
from pathlib import Path

import numpy as np
from matchms import Spectrum
from matchms.importing import load_spectra
from matchms.importing.load_from_json import as_spectrum

log = logging.getLogger(__name__)

# Whitespace and the commas between the records of a JSON list.
_BETWEEN_RECORDS = re.compile(r"[\s,]*")

# MoNA metadata names -> matchms keys. Compound metadata...
MONA_COMPOUND_FIELDS = {"molecular formula": "formula", "smiles": "smiles", "inchi": "inchi", "inchikey": "inchikey"}
# ...and spectrum metadata.
MONA_SPECTRUM_FIELDS = {
    "precursor m/z": "precursor_mz",
    "precursor type": "adduct",
    "ionization mode": "ionmode",
    "ms level": "ms_level",
    "retention time": "retention_time",
    "instrument": "instrument",
    "instrument type": "instrument_type",
    "collision energy": "collision_energy",
}


def read_library(path: str | Path) -> Iterator[Spectrum]:
    """The spectra of a library file (MSP, MGF or JSON), one at a time."""
    if Path(path).suffix.lower() == ".json":
        yield from read_json_library(path)
    else:
        yield from load_spectra(str(path))


# ---- JSON --------------------------------------------------------------------

def json_records(path: str | Path, chunk_chars: int = 1 << 22) -> Iterator[dict]:
    """
    The records of a JSON file that holds one list of records
    (``[{...}, {...}, ...]``), one at a time.

    The file is read in chunks of ``chunk_chars`` characters. Python's own
    JSON decoder (``raw_decode``) reads one record from the text at hand; when
    the record is cut by the end of the chunk, the decoder fails and the next
    chunk is appended.
    """
    decoder = json.JSONDecoder()
    with open(path, encoding="utf-8-sig") as handle:  # -sig: skip a byte order mark
        text = handle.read(chunk_chars)
        position = _BETWEEN_RECORDS.match(text, 0).end()
        if not text.startswith("[", position):
            raise ValueError(f"{Path(path).name}: not a JSON list of spectra")
        position += 1
        while True:
            position = _BETWEEN_RECORDS.match(text, position).end()
            if text.startswith("]", position):
                return  # end of the list
            try:
                record, position = decoder.raw_decode(text, position)
            except json.JSONDecodeError:
                more = handle.read(chunk_chars)
                if not more:
                    raise  # the file ends inside a record: truncated download?
                text, position = text[position:] + more, 0
                continue
            yield record


def read_json_library(path: str | Path) -> Iterator[Spectrum]:
    """The spectra of a GNPS, MoNA or MassBank JSON library."""
    n_unknown = 0
    for record in json_records(path):
        if "peaks_json" in record:
            spectrum = as_spectrum(record)  # GNPS
        elif "spectrum" in record and "compound" in record:
            spectrum = mona_spectrum(record)
        elif "PK$PEAK" in record:
            spectrum = massbank_spectrum(record)
        elif record.get("DEPRECATED"):
            continue  # MassBank keeps withdrawn records, without their peaks
        else:
            n_unknown += 1
            continue
        if spectrum is not None:
            yield spectrum
    if n_unknown:
        log.warning("%s: %d JSON records of an unknown format were skipped", Path(path).name, n_unknown)


def retention_minutes(value) -> float | None:
    """
    A retention time as written in MoNA / MassBank records ("5.68 min",
    "340.2 sec", "340 s", or a bare number, taken as minutes) in minutes.
    """
    if value is None:
        return None
    match = re.match(r"\s*(\d+(?:\.\d*)?)\s*([a-z]*)", str(value).lower())
    if match is None:
        return None
    number = float(match[1])
    return number / 60.0 if match[2].startswith("s") else number


def _spectrum(mz, intensities, metadata: dict) -> Spectrum | None:
    """A matchms spectrum from peaks and metadata (empty values dropped)."""
    if len(mz) == 0:
        return None
    order = np.argsort(mz)  # matchms needs sorted m/z
    metadata = {key: value for key, value in metadata.items() if value not in (None, "", [])}
    return Spectrum(mz=np.asarray(mz, dtype=float)[order], intensities=np.asarray(intensities, dtype=float)[order],
                    metadata=metadata)


def mona_spectrum(record: dict) -> Spectrum | None:
    """
    One record of a MoNA JSON export. Its tags and library name go in the
    comment, where the in-silico check looks for words such as "In-Silico"
    or "LipidBlast".
    """
    compound = (record.get("compound") or [{}])[0]
    names = compound.get("names") or [{}]
    metadata = {
        "spectrum_id": record.get("id"),
        "compound_name": names[0].get("name"),
        "inchi": compound.get("inchi"),
        "inchikey": compound.get("inchiKey"),
    }
    for item in compound.get("metaData") or []:
        key = MONA_COMPOUND_FIELDS.get(str(item.get("name")).lower())
        if key and not metadata.get(key):
            metadata[key] = item.get("value")
    for item in record.get("metaData") or []:
        key = MONA_SPECTRUM_FIELDS.get(str(item.get("name")).lower())
        if key:
            metadata[key] = item.get("value")
    if "retention_time" in metadata:
        metadata["retention_time"] = retention_minutes(metadata["retention_time"])
    tags = [str(tag.get("text", "")) for tag in record.get("tags") or []]
    library = (record.get("library") or {}).get("library")
    metadata["comment"] = "; ".join(filter(None, [*tags, library]))

    # Peaks as text: "184.0733:999 86.0964:12.5 ..."
    values = np.array(str(record.get("spectrum") or "").replace(":", " ").split(), dtype=float)
    peaks = values[: values.size // 2 * 2].reshape(-1, 2)
    return _spectrum(peaks[:, 0], peaks[:, 1], metadata)


def massbank_spectrum(record: dict) -> Spectrum | None:
    """One record of a MassBank JSON file (keys of the MassBank record format)."""
    focused_ion = record.get("MS$FOCUSED_ION") or {}
    chromatography = record.get("AC$CHROMATOGRAPHY") or {}
    metadata = {
        "spectrum_id": record.get("ACCESSION"),
        "compound_name": (record.get("CH$NAME") or [None])[0],
        "formula": record.get("CH$FORMULA"),
        "smiles": record.get("CH$SMILES"),
        "inchi": record.get("CH$IUPAC"),
        "inchikey": (record.get("CH$LINK") or {}).get("INCHIKEY"),
        "ionmode": record.get("AC$MASS_SPECTROMETRY_ION_MODE"),
        "ms_level": record.get("AC$MASS_SPECTROMETRY_MS_TYPE"),
        "instrument": record.get("AC$INSTRUMENT"),
        "instrument_type": record.get("AC$INSTRUMENT_TYPE"),
        "precursor_mz": focused_ion.get("PRECURSOR_M/Z"),
        "adduct": focused_ion.get("PRECURSOR_TYPE"),
        "retention_time": retention_minutes(chromatography.get("RETENTION_TIME")),
        "comment": "; ".join(record.get("COMMENT") or []),
    }
    # Peaks: [["m/z", "intensity", "relative intensity"], ...] as text.
    peaks = np.array(record.get("PK$PEAK") or [], dtype=float).reshape(-1, 3)
    return _spectrum(peaks[:, 0], peaks[:, 1], metadata)
