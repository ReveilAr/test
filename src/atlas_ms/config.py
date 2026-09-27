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

import re
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


class SpectrumQCSettings(Section):
    """
    Quality control of the MS2 spectra before similarity scoring. These
    filters are independent of the scores: never filter on a score itself.
    """

    precursor_window_da = param.Number(
        default=17.0, bounds=(0.0, 100.0),
        doc="Fragments within this distance (Da) of the precursor m/z are removed: "
            "the unfragmented precursor and its isotopes (GNPS default 17 Da).",
    )
    min_peaks = param.Integer(
        default=4, bounds=(1, None),
        doc="Spectra with fewer fragments (after cleaning) get no spectral edges. "
            "Lipid spectra are sparse (PC: m/z 184 and little else), so keep it low.",
    )
    max_blank_ratio = param.Number(
        default=0.3, bounds=(0.0, None),
        doc="Features whose mean intensity in blanks exceeds this fraction of their "
            "mean in samples get no spectral edges. Only applies when samples.tsv "
            "lists blanks (FBMN-STATS default 0.3).",
    )


class ScoringSettings(Section):
    """
    Pairwise spectral similarity. Scoring is the slow step, so it keeps a
    pool of candidate edges; the network is then built from the pool, and
    changing the network cutoffs does not re-run the scoring.
    """

    score = param.Selector(
        default="modified_cosine", objects=["modified_cosine", "ms2deepscore"],
        doc="modified_cosine: fragment matching, also across the precursor mass "
            "difference (GNPS). ms2deepscore: predicted structural similarity "
            "(pretrained deep-learning model).",
    )
    fragment_tolerance_da = param.Number(
        default=0.01, bounds=(0.0001, 1.0),
        doc="Modified cosine only: tolerance (Da) when matching fragments.",
    )
    ms2deepscore_model = param.String(
        default="",
        doc="MS2DeepScore only: path of the model file (.pt). Empty: the pretrained "
            "model, downloaded once to ~/.cache/atlas-ms/models/.",
    )
    candidates_per_spectrum = param.Integer(
        default=50, bounds=(1, None),
        doc="Best-scoring neighbours kept per spectrum as candidate edges. "
            "The network's top_k cannot be larger.",
    )
    min_candidate_score = param.Number(
        default=0.3, bounds=(0.0, 1.0),
        doc="Pairs scoring lower are never candidates. The network's min_score "
            "cannot be lower.",
    )


class NetworkSettings(Section):
    """
    Molecular network construction from the candidate edges (GNPS algorithm).
    Values are GNPS placeholders, to be tuned on the first real network.
    """

    min_score = param.Number(
        default=0.7, bounds=(0.0, 1.0),
        doc="Minimum similarity for an edge (placeholder: GNPS modified cosine value; "
            "MS2DeepScore usually needs a higher one).",
    )
    min_matched_peaks = param.Integer(
        default=6, bounds=(0, None),
        doc="Modified cosine only: minimum number of matched fragments for an edge. "
            "Lipid spectra are sparse; this GNPS value may be too strict for them.",
    )
    top_k = param.Integer(
        default=10, bounds=(1, None),
        doc="An edge is kept only if each node is among the other's top_k most "
            "similar neighbours.",
    )
    max_family_size = param.Integer(
        default=100, bounds=(0, None),
        doc="Families (connected components) larger than this lose their weakest "
            "edges until they split. 0: no limit.",
    )
    louvain_resolution = param.Number(
        default=1.0, bounds=(0.01, 10.0),
        doc="Louvain community detection inside families: higher values give "
            "smaller communities.",
    )


LIBRARY_KINDS = ("experimental", "in_silico")


class LibrarySearchSettings(Section):
    """
    Spectral library search (matchms): every feature's MS2 spectrum against
    the library spectra whose precursor m/z is within tolerance (identity
    search). Hits set the confidence level: 1 for reference standards with a
    matching retention time, 2a for other experimental spectra, 3 for
    in-silico (predicted) spectra.
    """

    libraries = param.List(
        default=[], item_type=dict,
        doc="Libraries as {name, path, kind, reference_standards, rt_unit, rt_tolerance_s} "
            "(edited in the app's library table). kind: experimental or in_silico.",
    )
    precursor_tolerance_ppm = param.Number(
        default=10.0, bounds=(0.1, 100.0),
        doc="A library spectrum is a candidate when its precursor m/z is within this tolerance.",
    )
    fragment_tolerance_da = param.Number(
        default=0.01, bounds=(0.0001, 1.0), doc="Tolerance (Da) when matching fragments.",
    )
    min_score = param.Number(
        default=0.7, bounds=(0.0, 1.0), doc="Minimum cosine score of a library hit.",
    )
    min_matched_peaks = param.Integer(
        default=3, bounds=(1, None),
        doc="Minimum number of matched fragments of a library hit. Lipid spectra are "
            "sparse, hence lower than the GNPS default (6).",
    )
    top_n = param.Integer(default=3, bounds=(1, None), doc="Hits kept per feature and library.")
    min_library_peaks = param.Integer(
        default=2, bounds=(1, None),
        doc="Library spectra with fewer fragments after cleaning are not used (lipid spectra are sparse).",
    )
    repair_annotations = param.Boolean(
        default=True,
        doc="Harmonization: derive missing SMILES / InChI / InChIKey / formula, repair annotations "
            "that contradict the precursor mass and remove those that cannot be repaired "
            "(spectra with a structure only). Slow for very large libraries (RDKit), but done once.",
    )
    remove_duplicates = param.Boolean(
        default=True, doc="Keep a spectrum found in several libraries once (from the most trusted library).",
    )

    def __init__(self, **params):
        super().__init__(**params)
        self.check_libraries()

    def check_libraries(self) -> None:
        """Complete entries, known kinds and units, unique file-name-safe names."""
        names = [library.get("name") for library in self.libraries]
        if len(set(names)) != len(names):
            raise ValueError(f"Library names must be unique: {names}")
        for library in self.libraries:
            missing = {"name", "path", "kind"} - set(library)
            if missing:
                raise ValueError(f"Library {library} is missing {sorted(missing)}")
            if not re.fullmatch(r"[A-Za-z0-9_-]+", str(library["name"])):
                raise ValueError(f"Library name '{library['name']}': use letters, digits, '_' or '-' only")
            if library["kind"] not in LIBRARY_KINDS:
                raise ValueError(f"Library '{library['name']}': kind must be one of {LIBRARY_KINDS}")
            if library.get("rt_unit", "min") not in ("min", "s"):
                raise ValueError(f"Library '{library['name']}': rt_unit must be 'min' or 's'")

    def library(self, name: str) -> dict:
        """
        One library's entry, with the optional keys filled in. Keys are
        sorted: Snakemake compares rule parameters as text, and project.yaml
        may store them in any order.
        """
        entry = next(library for library in self.libraries if library["name"] == name)
        defaults = {"reference_standards": False, "rt_unit": "min", "rt_tolerance_s": 10.0}
        return dict(sorted((defaults | entry).items()))


class LipidSettings(Section):
    """
    Rule-based lipid annotation: class-specific fragments and neutral losses
    (``atlas_ms/presets/lipid_rules.yaml``) plus the precursor m/z of every
    species of the class. Gives the species level (e.g. PC 34:1), or the
    chains (e.g. TG 16:0_18:1_18:1) when fragments reveal them.
    """

    enabled = param.Boolean(default=True)
    rules_file = param.String(
        default="",
        doc="Your own rules (YAML, same format as atlas_ms/presets/lipid_rules.yaml). "
            "Empty: the built-in rules.",
    )
    precursor_tolerance_ppm = param.Number(
        default=5.0, bounds=(0.1, 100.0), doc="Precursor m/z tolerance (ppm) for a lipid species.",
    )
    fragment_tolerance_da = param.Number(
        default=0.01, bounds=(0.0001, 1.0), doc="Tolerance (Da) for fragments and neutral losses.",
    )
    min_relative_intensity = param.Number(
        default=0.01, bounds=(0.0, 1.0),
        doc="Fragments weaker than this fraction of the base peak are ignored.",
    )
    mz_only_candidates = param.Boolean(
        default=True,
        doc="Also list the species matching the precursor m/z only, without fragment "
            "evidence (level 5 suggestions; the only ones for features without MS2).",
    )


class SiriusSettings(Section):
    """
    SIRIUS 6 (local REST service, installed by the pipeline in its own conda
    environment): molecular formula + ZODIAC, El Gordo lipids, CSI:FingerID
    structures with COSMIC confidence, CANOPUS classes. Off by default.
    """

    enabled = param.Boolean(
        default=False, label="Run SIRIUS",
        doc="Needs a free academic SIRIUS account for CSI:FingerID, CANOPUS and MSNovelist (below). "
            "Takes hours for thousands of features.",
    )
    profile = param.Selector(default="orbitrap", objects=["orbitrap", "qtof"], doc="SIRIUS instrument profile.")
    ms2_ppm = param.Number(default=5.0, bounds=(0.1, 50.0), label="MS2 mass accuracy (ppm)",
                           doc="MS2 mass accuracy (ppm) for formula identification.")
    formula_candidates = param.Integer(default=10, bounds=(1, None), doc="Formula candidates computed per feature.")
    zodiac = param.Boolean(default=True, label="Re-rank with ZODIAC",
                           doc="Re-rank formulas with ZODIAC (uses the whole data set).")
    min_formula_score = param.Number(
        default=0.9, bounds=(0.0, 1.0),
        doc="The top formula is level 4 when its ZODIAC score (SIRIUS score without ZODIAC) reaches this.",
    )
    structure_databases = param.List(
        default=["BIO"], item_type=str,
        doc="SIRIUS structure databases searched by CSI:FingerID (e.g. BIO, LIPIDMAPS, or your custom ones).",
    )
    expansive_search = param.Selector(
        default="APPROXIMATE", objects=["OFF", "EXACT", "APPROXIMATE"],
        doc="Search PubChem when the best database structure has a low COSMIC confidence "
            "(OFF, or by the exact / approximate confidence).",
    )
    canopus = param.Boolean(default=True, label="Run CANOPUS", doc="Predict compound classes (ClassyFire, NPClassifier).")
    msnovelist = param.Boolean(
        default=False, label="Run MSNovelist",
        doc="MSNovelist: generate structures de novo from the predicted fingerprint, for compounds "
            "missing from every database. Slow (a web service); candidates are level 3 at most.",
    )
    msnovelist_candidates = param.Integer(
        default=128, bounds=(1, None), label="MSNovelist candidates",
        doc="MSNovelist: structures generated per feature before ranking.",
    )
    candidates_kept = param.Integer(default=3, bounds=(1, None), doc="Candidates kept per feature and SIRIUS tool.")


class MS2QuerySettings(Section):
    """
    MS2Query analog search (own conda environment; its library and models,
    a few GB, are downloaded once per machine). Off by default.
    """

    enabled = param.Boolean(default=False)
    models_dir = param.String(
        default="", doc="Folder with MS2Query's positive-mode library and models. Empty: downloaded "
                        "once to ~/.cache/atlas-ms/ms2query/positive/.",
    )
    top_n = param.Integer(default=3, bounds=(1, None), doc="Analogs kept per feature.")
    min_score = param.Number(
        default=0.7, bounds=(0.0, 1.0), doc="Minimum MS2Query score (the authors report 0.7 as reliable).",
    )
    precursor_tolerance_ppm = param.Number(
        default=10.0, bounds=(0.1, 100.0), doc="An analog within this m/z of the feature is marked as an exact match.",
    )


class HarmonizationSettings(Section):
    """
    Combination of all annotation sources: confidence levels (Schymanski),
    flags, best annotation per feature, class consensus per family.
    """

    source_priority = param.List(
        default=["library", "lipid_rules", "sirius:elgordo", "ms2query", "sirius:csi", "sirius:msnovelist",
                 "sirius:canopus", "sirius"],
        item_type=str,
        doc="Between candidates of the same level, the best comes from the first source in this list "
            "(a name like 'library' covers every library; 'sirius:csi' one SIRIUS tool).",
    )
    rt_model_min_points = param.Integer(
        default=5, bounds=(3, None),
        doc="A lipid class needs this many confident annotations (level 3 or better) "
            "before its retention-time trend (RT vs carbons and double bonds) is used.",
    )
    rt_outlier_min = param.Number(
        default=0.5, bounds=(0.0, None),
        doc="Annotations further than this (minutes) from their class's RT trend are flagged.",
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
    "spectrum_qc": SpectrumQCSettings,
    "scoring": ScoringSettings,
    "network": NetworkSettings,
    "library_search": LibrarySearchSettings,
    "lipids": LipidSettings,
    "sirius": SiriusSettings,
    "ms2query": MS2QuerySettings,
    "harmonization": HarmonizationSettings,
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
    spectrum_qc: SpectrumQCSettings
    scoring: ScoringSettings
    network: NetworkSettings
    library_search: LibrarySearchSettings
    lipids: LipidSettings
    sirius: SiriusSettings
    ms2query: MS2QuerySettings
    harmonization: HarmonizationSettings

    def __init__(self, **sections: Section):
        unknown = set(sections) - set(SECTIONS)
        if unknown:
            raise ValueError(f"Unknown section(s): {sorted(unknown)}")
        for key, section_class in SECTIONS.items():
            setattr(self, key, sections.get(key) or section_class())

    def check(self) -> None:
        """Rules that involve several sections (single values are checked by param)."""
        if self.network.min_score < self.scoring.min_candidate_score:
            raise ValueError(
                f"network.min_score ({self.network.min_score}) is below "
                f"scoring.min_candidate_score ({self.scoring.min_candidate_score}): "
                "lower the candidate score, or raise the network one"
            )
        if self.network.top_k > self.scoring.candidates_per_spectrum:
            raise ValueError(
                f"network.top_k ({self.network.top_k}) is larger than "
                f"scoring.candidates_per_spectrum ({self.scoring.candidates_per_spectrum})"
            )

    # ---- dict / YAML round trip ----

    def to_dict(self) -> dict[str, dict]:
        return {key: getattr(self, key).to_dict() for key in SECTIONS}

    @classmethod
    def from_dict(cls, data: dict | None) -> "ProjectConfig":
        data = dict(data or {})
        unknown = set(data) - set(SECTIONS)
        if unknown:
            raise ValueError(f"Unknown section(s) in project configuration: {sorted(unknown)}")
        config = cls(**{key: SECTIONS[key].from_dict(data.get(key)) for key in SECTIONS})
        config.check()
        return config

    @classmethod
    def load(cls, path: str | Path) -> "ProjectConfig":
        return cls.from_dict(yaml.safe_load(Path(path).read_text()))

    def save(self, path: str | Path) -> None:
        self.check()
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
