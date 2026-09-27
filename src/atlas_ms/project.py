"""
Project folders and the sample table.

A project is a folder that holds everything about one study:

    my_project/
    ├── project.yaml   all parameters (see atlas_ms.config.ProjectConfig)
    ├── samples.tsv    one row per raw file (see below)
    ├── work/          intermediate files, managed by Snakemake
    ├── results/       final tables and exports, read by the app
    └── logs/          one log file per processing step and sample

``samples.tsv`` columns:

    sample       short unique name, used in every file name (letters, digits, _ and -)
    file         absolute path of the raw file (.raw from Thermo, or .mzML)
    sample_type  "sample", "blank" or "standard"
    ATTRIBUTE_*  any number of free metadata columns (groups, time points...),
                 named as in the GNPS metadata format

Raw files are referenced by path, never copied into the project.
"""

import re
from pathlib import Path

import pandas as pd

from atlas_ms.config import ProjectConfig

SAMPLE_TYPES = ("sample", "blank", "standard")

# Raw file extensions the pipeline can read, and how each one is prepared.
# Thermo .raw files are converted with ThermoRawFileParser; mzML files are
# used as they are (they must already be centroided).
RAW_FORMATS = {".raw": "thermo", ".mzml": "mzml"}


def raw_format(path: str | Path) -> str:
    """Return "thermo" or "mzml" for a raw file, or raise for unsupported files."""
    suffix = Path(path).suffix.lower()
    if suffix not in RAW_FORMATS:
        raise ValueError(f"Unsupported raw file '{path}': expected one of {sorted(RAW_FORMATS)}")
    return RAW_FORMATS[suffix]


def sample_name(path: str | Path) -> str:
    """File name without extension, reduced to characters safe in file names and Snakemake."""
    return re.sub(r"[^A-Za-z0-9_-]+", "_", Path(path).stem).strip("_") or "sample"


def make_sample_table(files: list[str | Path], taken=frozenset()) -> pd.DataFrame:
    """
    Build a new sample table (all rows of type "sample") from raw file paths.
    Sample names are unique, also among ``taken`` (names already in use).
    """
    rows, used = [], set(taken)
    for file in files:
        raw_format(file)  # fail early on unsupported files
        name = sample_name(file)
        # Two files can share a stem (e.g. same name in different folders).
        unique, n = name, 2
        while unique in used:
            unique, n = f"{name}_{n}", n + 1
        used.add(unique)
        rows.append({"sample": unique, "file": str(Path(file).resolve()), "sample_type": "sample"})
    return pd.DataFrame(rows, columns=["sample", "file", "sample_type"])


def add_samples(samples: pd.DataFrame, files: list[str | Path]) -> pd.DataFrame:
    """
    The sample table with a new row (type "sample", empty metadata) for each
    file not listed yet. Adding samples to a processed project re-runs the
    steps that combine all runs (alignment, linking...) on the next run.
    """
    listed = {str(Path(f).resolve()) for f in samples["file"]}
    new = make_sample_table([f for f in files if str(Path(f).resolve()) not in listed], taken=set(samples["sample"]))
    return pd.concat([samples, new], ignore_index=True).fillna("")


def validate_samples(samples: pd.DataFrame) -> pd.DataFrame:
    """
    Check a sample table and return it indexed by sample name.

    Raises ValueError with a readable message on the first problem found.
    """
    for column in ("sample", "file", "sample_type"):
        if column not in samples.columns:
            raise ValueError(f"samples.tsv has no '{column}' column")
    if samples.empty:
        raise ValueError("samples.tsv lists no raw files")
    duplicated = samples["sample"][samples["sample"].duplicated()].tolist()
    if duplicated:
        raise ValueError(f"Duplicated sample names: {duplicated}")
    bad_names = [s for s in samples["sample"] if s != sample_name(s)]
    if bad_names:
        raise ValueError(f"Sample names may only contain letters, digits, '_' and '-': {bad_names}")
    bad_types = sorted(set(samples["sample_type"]) - set(SAMPLE_TYPES))
    if bad_types:
        raise ValueError(f"Unknown sample_type {bad_types}; use one of {SAMPLE_TYPES}")
    other = [c for c in samples.columns if c not in ("sample", "file", "sample_type")]
    not_attribute = [c for c in other if not c.startswith("ATTRIBUTE_")]
    if not_attribute:
        raise ValueError(f"Metadata columns must start with 'ATTRIBUTE_': {not_attribute}")
    for file in samples["file"]:
        raw_format(file)
    return samples.set_index("sample", drop=False)


def read_samples(path: str | Path) -> pd.DataFrame:
    """Read and validate ``samples.tsv``; the result is indexed by sample name."""
    samples = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return validate_samples(samples)


def write_samples(samples: pd.DataFrame, path: str | Path) -> None:
    validate_samples(samples.reset_index(drop=True))
    samples.to_csv(path, sep="\t", index=False)


class Project:
    """Paths and file access for one project folder."""

    def __init__(self, root: str | Path):
        self.root = Path(root).resolve()

    @property
    def config_path(self) -> Path:
        return self.root / "project.yaml"

    @property
    def samples_path(self) -> Path:
        return self.root / "samples.tsv"

    @property
    def results_dir(self) -> Path:
        return self.root / "results"

    def load_config(self) -> ProjectConfig:
        return ProjectConfig.load(self.config_path)

    def save_config(self, config: ProjectConfig) -> None:
        config.save(self.config_path)

    def load_samples(self) -> pd.DataFrame:
        return read_samples(self.samples_path)

    def save_samples(self, samples: pd.DataFrame) -> None:
        write_samples(samples, self.samples_path)

    @classmethod
    def create(
        cls,
        root: str | Path,
        files: list[str | Path],
        instrument: str = "orbitrap",
        adducts: str = "positive_lipids",
    ) -> "Project":
        """
        Create a new project folder with default parameters (from the given
        presets) and a sample table listing ``files``. Refuses to overwrite an
        existing project.
        """
        project = cls(root)
        if project.config_path.exists():
            raise FileExistsError(f"{project.root} already contains a project")
        project.root.mkdir(parents=True, exist_ok=True)
        config = ProjectConfig()
        config.apply_instrument_preset(instrument)
        config.apply_adduct_preset(adducts)
        project.save_config(config)
        project.save_samples(make_sample_table(files))
        return project
