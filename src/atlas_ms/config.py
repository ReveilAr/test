"""
Project parameters.

Every tunable value of the pipeline is declared here, once, as an attribute
of a `param.Parameterized` "section" (one section per processing step):

* the Panel app renders each section as a group of widgets (``pn.Param``);
* ``project.yaml`` is the serialised form of all sections (``ProjectConfig``);
* the Snakefile rebuilds the sections from ``project.yaml``, so wrong types or
  out-of-range values are rejected before any rule runs;
* each rule receives the plain-dict form of the sections it uses, so
  Snakemake re-runs a rule when (and only when) one of its values changes.

Presets (``atlas_ms/presets/*.yaml``) are named sets of values that can be
copied into a project in one go: instrument presets set tolerances across
several sections, adduct presets set the polarity and the adduct list.
"""

from importlib import resources
from pathlib import Path
from typing import Any

import param
import yaml


# --------------------------------------------------------------------------
# Presets
# --------------------------------------------------------------------------

def load_presets(kind: str) -> dict[str, dict]:
    """Read ``atlas_ms/presets/<kind>.yaml`` (``kind`` is "instruments" or "adducts")."""
    text = (resources.files("atlas_ms") / "presets" / f"{kind}.yaml").read_text()
    return yaml.safe_load(text)


def _default_adducts() -> list[dict]:
    """Adduct list of the default preset (lipidomics, positive mode with ammonium)."""
    return load_presets("adducts")["positive_lipids"]["adducts"]


# --------------------------------------------------------------------------
# Sections
# --------------------------------------------------------------------------

class Section(param.Parameterized):
    """
    Base class of all parameter sections.

    Adds conversion to and from plain dicts, which is what ``project.yaml``
    stores and what Snakemake passes to the rule scripts.
    """

    def to_dict(self) -> dict[str, Any]:
        values = self.param.values()
        values.pop("name")  # param's internal object name, not a setting
        return values

    @classmethod
    def from_dict(cls, values: dict[str, Any] | None) -> "Section":
        values = dict(values or {})
        # Reject unknown keys so that a typo in project.yaml is an error,
        # not a silently ignored setting. Missing keys keep their defaults.
        unknown = set(values) - set(cls.param) - {"name"}
        if unknown:
            raise ValueError(f"Unknown parameter(s) for '{cls.__name__}': {sorted(unknown)}")
        return cls(**values)


class PresetInfo(Section):
    """
    Names of the presets the parameters were initialised from. For display
    only: no processing step reads them, so renaming never triggers a re-run.
    """

    instrument = param.String(default="orbitrap")
    adducts = param.String(default="positive_lipids")


class InstrumentSettings(Section):
    """
    Values that depend on the mass spectrometer and the chromatography,
    used by feature finding (and by gap filling for the peak width).
    """

    mass_error_ppm = param.Number(
        default=10.0, bounds=(0.1, 100.0),
        doc="Allowed m/z deviation (ppm) of consecutive centroids of one mass trace.",
    )
    noise_threshold = param.Number(
        default=1e4, bounds=(0.0, None),
        doc="MS1 centroids below this intensity are treated as noise.",
    )
    chrom_fwhm_s = param.Number(
        default=5.0, bounds=(0.1, 120.0),
        doc="Typical chromatographic peak width at half height (seconds).",
    )
    min_trace_length_s = param.Number(
        default=5.0, bounds=(0.0, None),
        doc="Mass traces shorter than this (seconds) are discarded.",
    )


class AdductSettings(Section):
    """Ion polarity and the adducts used to group ions of the same molecule."""

    polarity = param.Selector(default="positive", objects=["positive", "negative"])
    adducts = param.List(
        default=_default_adducts(), item_type=dict,
        doc="Adducts as {name, openms, probability} (see atlas_ms/presets/adducts.yaml).",
    )

    def __init__(self, **params):
        super().__init__(**params)
        self._check_adducts()

    def _check_adducts(self) -> None:
        """Each adduct must be complete, and its charge must match the polarity."""
        wrong_sign = "-" if self.polarity == "positive" else "+"
        for adduct in self.adducts:
            missing = {"name", "openms", "probability"} - set(adduct)
            if missing:
                raise ValueError(f"Adduct {adduct} is missing {sorted(missing)}")
            charge = adduct["openms"].rsplit(":", 1)[-1]
            if charge.startswith(wrong_sign):
                raise ValueError(f"Adduct {adduct['name']} does not match {self.polarity} polarity")
            if adduct["probability"] <= 0:
                raise ValueError(f"Adduct {adduct['name']} needs a probability > 0")

    def openms_adducts(self) -> list[bytes]:
        """
        Adducts in the "<formula>:<charge>:<probability>" format of
        MetaboliteAdductDecharger. OpenMS requires the probabilities of the
        charged adducts to sum to 1, so they are rescaled here; neutral gains
        and losses (charge "0") keep their own probability.
        """
        def is_neutral(adduct: dict) -> bool:
            return adduct["openms"].endswith(":0")

        total = sum(a["probability"] for a in self.adducts if not is_neutral(a))
        result = []
        for adduct in self.adducts:
            p = adduct["probability"] if is_neutral(adduct) else adduct["probability"] / total
            result.append(f"{adduct['openms']}:{p:.4g}".encode())
        return result


class FeatureFindingSettings(Section):
    """Untargeted feature detection in each run (pyOpenMS FeatureFinderMetabo)."""

    remove_single_traces = param.Boolean(
        default=True,
        doc="Drop features made of a single mass trace (no isotope peak detected). "
            "Removes a lot of noise, but also compounds so weak that their M+1 "
            "peak is below the noise threshold.",
    )
    isotope_filtering_model = param.Selector(
        default="none",
        objects=["none", "metabolites (2% RMS)", "metabolites (5% RMS)"],
        doc="Check isotope patterns against a model of average metabolites. "
            "'none' (UmetaFlow default) keeps unusual patterns, e.g. of Cl/Br compounds.",
    )
    precursor_peak_tolerance_ppm = param.Number(
        default=100.0, bounds=(0.0, 1000.0),
        doc="Before feature finding, each MS2 precursor m/z is moved to the most "
            "intense MS1 peak within this window (corrects the isolation-window centre).",
    )
    precursor_feature_tolerance_ppm = param.Number(
        default=5.0, bounds=(0.0, 100.0),
        doc="After feature finding, each MS2 precursor m/z is moved to the m/z of "
            "the feature it belongs to, if within this tolerance.",
    )


class AlignmentSettings(Section):
    """Retention-time alignment of the runs (MapAlignerPoseClustering)."""

    mz_max_ppm = param.Number(
        default=10.0, bounds=(0.1, 100.0),
        doc="Features further apart than this in m/z are never paired between runs.",
    )
    rt_max_difference_s = param.Number(
        default=100.0, bounds=(1.0, None),
        doc="Largest retention-time difference (s) between paired features of two runs.",
    )


class LinkingSettings(Section):
    """Grouping of corresponding features across runs (FeatureLinkerUnlabeledKD)."""

    mz_tol_ppm = param.Number(default=8.0, bounds=(0.1, 100.0), doc="m/z tolerance (ppm).")
    rt_tol_s = param.Number(default=30.0, bounds=(0.1, None), doc="Retention-time tolerance (s), after alignment.")


class GapFillingSettings(Section):
    """
    Targeted re-extraction of features missing in some runs
    (FeatureFinderMetaboIdent), so every sample gets a value.
    """

    enabled = param.Boolean(default=True)
    mz_window_ppm = param.Number(default=10.0, bounds=(0.1, 100.0), doc="Extraction window (ppm).")
    rt_window_s = param.Number(default=30.0, bounds=(1.0, None), doc="Extraction window (s) around the expected retention time.")


class ExportSettings(Section):
    """Final feature table and GNPS/FBMN export."""

    min_detection_fraction = param.Number(
        default=0.0, bounds=(0.0, 1.0),
        doc="Keep a feature only if found in at least this fraction of the runs "
            "(after gap filling). 0 keeps everything.",
    )
    ms2_spectrum = param.Selector(
        default="most_intense", objects=["most_intense", "merged_spectra"],
        doc="MS2 spectrum exported per feature: the most intense one, or a merge "
            "of all similar MS2 spectra of the feature across runs.",
    )


# --------------------------------------------------------------------------
# All sections of a project
# --------------------------------------------------------------------------

# Section name in project.yaml -> section class, in processing order.
SECTIONS: dict[str, type[Section]] = {
    "presets": PresetInfo,
    "instrument": InstrumentSettings,
    "adducts": AdductSettings,
    "feature_finding": FeatureFindingSettings,
    "alignment": AlignmentSettings,
    "linking": LinkingSettings,
    "gap_filling": GapFillingSettings,
    "export": ExportSettings,
}


class ProjectConfig:
    """
    All parameters of a project: one section object per processing step,
    available as attributes (``config.instrument.noise_threshold``).
    """

    presets: PresetInfo
    instrument: InstrumentSettings
    adducts: AdductSettings
    feature_finding: FeatureFindingSettings
    alignment: AlignmentSettings
    linking: LinkingSettings
    gap_filling: GapFillingSettings
    export: ExportSettings

    def __init__(self, **sections: Section):
        unknown = set(sections) - set(SECTIONS)
        if unknown:
            raise ValueError(f"Unknown section(s): {sorted(unknown)}")
        for key, section_class in SECTIONS.items():
            setattr(self, key, sections.get(key) or section_class())

    # ---- dict / YAML round trip ----

    def to_dict(self) -> dict[str, dict]:
        return {key: getattr(self, key).to_dict() for key in SECTIONS}

    @classmethod
    def from_dict(cls, data: dict | None) -> "ProjectConfig":
        data = dict(data or {})
        unknown = set(data) - set(SECTIONS)
        if unknown:
            raise ValueError(f"Unknown section(s) in project configuration: {sorted(unknown)}")
        return cls(**{key: SECTIONS[key].from_dict(data.get(key)) for key in SECTIONS})

    @classmethod
    def load(cls, path: str | Path) -> "ProjectConfig":
        return cls.from_dict(yaml.safe_load(Path(path).read_text()))

    def save(self, path: str | Path) -> None:
        header = "# ATLAS-MS project parameters (edit here or in the app).\n"
        Path(path).write_text(header + yaml.safe_dump(self.to_dict(), sort_keys=False))

    # ---- presets ----

    def apply_instrument_preset(self, name: str) -> None:
        """Copy the values of an instrument preset into the matching sections."""
        presets = load_presets("instruments")
        if name not in presets:
            raise ValueError(f"Unknown instrument preset '{name}'. Available: {sorted(presets)}")
        for key, values in presets[name].items():
            if key == "description":
                continue
            getattr(self, key).param.update(**values)
        self.presets.instrument = name

    def apply_adduct_preset(self, name: str) -> None:
        """Replace polarity and adduct list with those of an adduct preset."""
        presets = load_presets("adducts")
        if name not in presets:
            raise ValueError(f"Unknown adduct preset '{name}'. Available: {sorted(presets)}")
        preset = presets[name]
        self.adducts = AdductSettings(polarity=preset["polarity"], adducts=preset["adducts"])
        self.presets.adducts = name
