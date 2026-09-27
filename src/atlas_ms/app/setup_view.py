"""
The Setup tab: sample metadata and processing parameters of the open project.

* Sample table: edit sample types (sample / blank / standard) and free
  metadata columns (saved as ATTRIBUTE_<name>, the GNPS convention); add
  samples by dropping raw files, remove selected ones.
* Presets: copy an instrument or adduct preset into the parameters.
* Adducts: polarity and the editable adduct list.
* Spectral libraries: the MSP / MGF / JSON files to search (typed paths or
  dropped files), with their kind (experimental or in-silico) and whether
  they are reference standards measured on the same method (the only way to
  reach level 1).
* Models: drop an MS2DeepScore model or MS2Query's files instead of
  downloading them.
* Parameters: every section of ``atlas_ms.config``, as input boxes; the
  SIRIUS card also holds the SIRIUS account.

Dropped files are copies (see ``uploads.py``).

Nothing is written until "Save" (or Run, which saves first). The files are
validated by the same code as the pipeline, so an invalid value is reported
here rather than failing a run.
"""

import re
from pathlib import Path

import pandas as pd
import panel as pn
import param

from atlas_ms.app.uploads import drop_area
from atlas_ms.config import LIBRARY_KINDS, SECTIONS, AdductSettings, load_presets
from atlas_ms.credentials import ACCOUNT_FILE, forget_sirius_account, load_sirius_account, save_sirius_account
from atlas_ms.network.model_files import DEFAULT_MODEL_PATH
from atlas_ms.project import SAMPLE_TYPES, Project, add_samples

# Sections shown as plain parameter boxes (presets and adducts have their own widgets).
PARAMETER_SECTIONS = [name for name in SECTIONS if name not in ("presets", "adducts")]
# Card titles that differ from the capitalised section name.
TITLES = {"spectrum_qc": "Spectrum QC", "lipids": "Lipid rules", "sirius": "SIRIUS", "ms2query": "MS2Query"}
HIDDEN = {"name", "libraries"}  # parameters with their own editor (libraries: the library table)
LIBRARY_COLUMNS = ["name", "path", "kind", "reference_standards", "rt_unit", "rt_tolerance_s"]
MS2QUERY_DROP_DIR = Path.home() / ".cache" / "atlas-ms" / "ms2query" / "dropped"
# The SIRIUS card, by topic (every parameter of SiriusSettings appears once).
SIRIUS_GROUPS = {
    "Formulas (SIRIUS, ZODIAC, El Gordo)": ["profile", "ms2_ppm", "formula_candidates", "zodiac", "min_formula_score"],
    "Structures (CSI:FingerID)": ["structure_databases", "expansive_search", "candidates_kept"],
    "De novo structures (MSNovelist)": ["msnovelist", "msnovelist_candidates"],
    "Compound classes (CANOPUS)": ["canopus"],
}

# Widget used for each kind of parameter: input boxes rather than sliders,
# so exact values (e.g. 10000.0) are easy to type.
WIDGETS = {
    param.Number: pn.widgets.FloatInput,
    param.Integer: pn.widgets.IntInput,
    param.Boolean: pn.widgets.Checkbox,
    param.Selector: pn.widgets.Select,
    param.String: pn.widgets.TextInput,
}


def widgets_for(section: param.Parameterized) -> dict:
    """Widget type for every parameter of a section (Integer before Number: it is a subclass)."""
    chosen = {}
    for name, parameter in section.param.objects("existing").items():
        if name == "name":
            continue
        for kind in (param.Integer, param.Number, param.Boolean, param.Selector, param.String):
            if isinstance(parameter, kind):
                chosen[name] = WIDGETS[kind]
                break
    return chosen


class SetupTab:
    """Widgets and actions of the Setup tab."""

    def __init__(self, project: Project, account_file: Path = ACCOUNT_FILE):
        self.project = project
        self.config = project.load_config()
        self.account_file = account_file  # the SIRIUS account (outside the project)
        self.message = pn.pane.Alert("", alert_type="light", visible=False, sizing_mode="stretch_width")

        # ---- sample table ----
        self.samples = pn.widgets.Tabulator(
            project.load_samples().reset_index(drop=True), selectable=True,
            show_index=False, sizing_mode="stretch_width", layout="fit_data_stretch",
            # Sample names and files are fixed; types come from a list; metadata is free text.
            editors={"sample": None, "file": None,
                     "sample_type": {"type": "list", "values": list(SAMPLE_TYPES)}},
        )
        self.new_column = pn.widgets.TextInput(placeholder="new metadata column, e.g. species", width=260)
        add_column = pn.widgets.Button(label="Add column", width=110)
        add_column.on_click(self._add_column)
        remove_samples = pn.widgets.Button(label="Remove selected samples", width=190)
        remove_samples.on_click(self._remove_samples)
        self.sample_drop = drop_area(str(project.root / "raw"), [".raw", ".mzml"],
                                     "Drop .raw / .mzML files here to add samples (copied into the project's raw/ folder)")
        self.sample_drop.param.watch(self._samples_dropped, "saved")

        # ---- presets ----
        self.instrument_preset = pn.widgets.Select(
            label="Instrument preset", options=list(load_presets("instruments")), value=self.config.presets.instrument)
        self.adduct_preset = pn.widgets.Select(
            label="Adduct preset", options=list(load_presets("adducts")), value=self.config.presets.adducts)
        apply_instrument = pn.widgets.Button(label="Apply", width=80)
        apply_adducts = pn.widgets.Button(label="Apply", width=80)
        apply_instrument.on_click(self._apply_instrument_preset)
        apply_adducts.on_click(self._apply_adduct_preset)

        # ---- adducts ----
        self.polarity = pn.widgets.Select(label="Polarity", options=["positive", "negative"],
                                          value=self.config.adducts.polarity, width=200)
        self.adducts = pn.widgets.Tabulator(
            pd.DataFrame(self.config.adducts.adducts, columns=["name", "openms", "probability"]),
            show_index=False, sizing_mode="stretch_width", selectable=True,
        )
        add_adduct = pn.widgets.Button(label="Add adduct", width=110)
        remove_adducts = pn.widgets.Button(label="Remove selected", width=140)
        add_adduct.on_click(self._add_adduct)
        remove_adducts.on_click(self._remove_adducts)

        # ---- spectral libraries ----
        self.libraries = pn.widgets.Tabulator(
            pd.DataFrame([self.config.library_search.library(lib["name"]) for lib in self.config.library_search.libraries],
                         columns=LIBRARY_COLUMNS),
            show_index=False, sizing_mode="stretch_width", selectable=True,
            titles={"reference_standards": "reference standards", "rt_unit": "library RT unit",
                    "rt_tolerance_s": "RT tolerance (s)"},
            editors={"kind": {"type": "list", "values": list(LIBRARY_KINDS)},
                     "reference_standards": {"type": "tickCross"},
                     "rt_unit": {"type": "list", "values": ["min", "s"]}},
            formatters={"reference_standards": {"type": "tickCross"}},
        )
        add_library = pn.widgets.Button(label="Add library", width=110)
        remove_libraries = pn.widgets.Button(label="Remove selected", width=140)
        add_library.on_click(self._add_library)
        remove_libraries.on_click(self._remove_libraries)
        self.library_drop = drop_area(str(project.root / "libraries"), [".msp", ".mgf", ".json"],
                                      "Drop MSP / MGF / JSON libraries here (copied into the project's libraries/ folder)")
        self.library_drop.param.watch(self._libraries_dropped, "saved")
        self.library_summary = pn.pane.DataFrame(self._read_library_summary(), index=False, sizing_mode="stretch_width")

        # ---- models ----
        self.ms2deepscore_drop = drop_area(str(DEFAULT_MODEL_PATH.parent), [".pt"],
                                           "Drop an MS2DeepScore model (.pt)", multiple=False)
        self.ms2deepscore_drop.param.watch(self._model_dropped, "saved")
        self.ms2query_drop = drop_area(str(MS2QUERY_DROP_DIR), [],
                                       "Drop MS2Query's library and model files (all of them)")
        self.ms2query_drop.param.watch(self._ms2query_dropped, "saved")
        for dropper in (self.sample_drop, self.library_drop, self.ms2deepscore_drop, self.ms2query_drop):
            dropper.param.watch(self._files_rejected, "rejected")

        # ---- SIRIUS account ----
        account = load_sirius_account(self.account_file) or {}
        self.sirius_user = pn.widgets.TextInput(label="SIRIUS account (e-mail)", value=account.get("username", ""),
                                                width=300)
        self.sirius_password = pn.widgets.PasswordInput(
            label="Password", placeholder="saved" if account else "", width=300)
        self.sirius_terms = pn.widgets.Checkbox(
            label="I accept the SIRIUS (Bright Giant) terms of service: academic, non-commercial use",
            value=account.get("accept_terms", False))
        save_login = pn.widgets.Button(label="Save login", width=110)
        forget_login = pn.widgets.Button(label="Forget login", width=110)
        save_login.on_click(lambda event: self.save_sirius_login())
        forget_login.on_click(lambda event: self.forget_sirius_login())
        self.sirius_status = pn.pane.Markdown(self._sirius_status_text())
        self.sirius_account = pn.Column(
            pn.pane.Markdown("**Account.** CSI:FingerID, CANOPUS and MSNovelist need a (free, academic) SIRIUS "
                             "account. It is stored for you only, outside the projects, and used when SIRIUS is "
                             "not logged in yet (it is checked at the next SIRIUS run)."),
            self.sirius_user, self.sirius_password, self.sirius_terms, pn.Row(save_login, forget_login),
            self.sirius_status,
        )

        save = pn.widgets.Button(label="Save", color="primary", width=120)
        save.on_click(lambda event: self.save())

        # One card per parameter section, bound to the config objects: editing
        # a box changes the section directly, applying a preset updates the boxes.
        self._cards = pn.Column(*self._parameter_cards(), sizing_mode="stretch_width")
        self.layout = pn.Column(
            self.message,
            pn.pane.Markdown("### Samples"),
            self.samples,
            pn.Row(self.new_column, add_column, remove_samples),
            self.sample_drop,
            pn.pane.Markdown("### Presets"),
            pn.Row(self.instrument_preset, apply_instrument, self.adduct_preset, apply_adducts),
            pn.pane.Markdown("### Adducts"),
            self.polarity,
            self.adducts,
            pn.Row(add_adduct, remove_adducts),
            pn.pane.Markdown("### Spectral libraries\n"
                             "MSP, MGF or GNPS JSON files (full paths). Only libraries of reference standards "
                             "measured on your method, with retention times, can give level 1."),
            self.libraries,
            pn.Row(add_library, remove_libraries),
            self.library_drop,
            pn.pane.Markdown("Last harmonization (spectra read, cleaned, marked in silico, duplicates removed, "
                             "searched; details in `results/annotations/library_cleaning.tsv`):"),
            self.library_summary,
            pn.pane.Markdown("### Models\nInstead of the automatic downloads (e.g. without internet access)."),
            pn.Row(self.ms2deepscore_drop, self.ms2query_drop, sizing_mode="stretch_width"),
            pn.pane.Markdown("### Parameters"),
            self._cards,
            save,
            sizing_mode="stretch_width",
        )

    def _parameter_cards(self) -> list:
        cards = []
        for name in PARAMETER_SECTIONS:
            section = getattr(self.config, name)
            if name == "sirius":  # on/off, the account, then the settings of each tool
                content = [pn.Param(section, parameters=["enabled"], widgets=widgets_for(section), show_name=False),
                           self.sirius_account]
                for title, parameters in SIRIUS_GROUPS.items():
                    content += [pn.pane.Markdown(f"**{title}**"),
                                pn.Param(section, parameters=parameters, widgets=widgets_for(section), show_name=False)]
            else:
                parameters = [p for p in section.param if p not in HIDDEN]
                content = [pn.Param(section, parameters=parameters, widgets=widgets_for(section), show_name=False)]
            cards.append(pn.Card(*content, title=TITLES.get(name, name.replace("_", " ").capitalize()),
                                 collapsed=True, sizing_mode="stretch_width"))
        return cards

    # ---- actions ----

    def _say(self, text: str, kind: str = "success") -> None:
        self.message.object, self.message.alert_type, self.message.visible = text, kind, True

    def _add_column(self, event) -> None:
        name = re.sub(r"[^A-Za-z0-9_]+", "_", self.new_column.value.strip()).strip("_")
        if not name:
            return
        column = name if name.startswith("ATTRIBUTE_") else f"ATTRIBUTE_{name}"
        self.samples.value = self.samples.value.assign(**{column: ""})
        self.new_column.value = ""

    def _apply_instrument_preset(self, event) -> None:
        self.config.apply_instrument_preset(self.instrument_preset.value)
        self._say(f"Instrument preset '{self.instrument_preset.value}' applied (not saved yet).", "info")

    def _apply_adduct_preset(self, event) -> None:
        self.config.apply_adduct_preset(self.adduct_preset.value)
        self.polarity.value = self.config.adducts.polarity
        self.adducts.value = pd.DataFrame(self.config.adducts.adducts, columns=["name", "openms", "probability"])
        self._say(f"Adduct preset '{self.adduct_preset.value}' applied (not saved yet).", "info")

    def _add_adduct(self, event) -> None:
        row = pd.DataFrame([{"name": "[M+?]+", "openms": "H:+", "probability": 0.1}])
        self.adducts.value = pd.concat([self.adducts.value, row], ignore_index=True)

    def _remove_adducts(self, event) -> None:
        self.adducts.value = self.adducts.value.drop(index=self.adducts.selection).reset_index(drop=True)
        self.adducts.selection = []

    # ---- dropped files ----

    @staticmethod
    def _new_files(event) -> list[str]:
        return event.new[len(event.old):]  # `saved` only grows

    def _samples_dropped(self, event) -> None:
        files = self._new_files(event)
        self.samples.value = add_samples(self.samples.value, files)
        self._say(f"{len(files)} sample(s) added (not saved yet): {', '.join(Path(f).name for f in files)}.", "info")

    def _remove_samples(self, event) -> None:
        self.samples.value = self.samples.value.drop(index=self.samples.selection).reset_index(drop=True)
        self.samples.selection = []

    def _libraries_dropped(self, event) -> None:
        table = self.libraries.value
        for path in self._new_files(event):
            name = re.sub(r"[^A-Za-z0-9_-]+", "_", Path(path).stem).strip("_") or "library"
            unique, n = name, 2
            while unique in set(table["name"]):
                unique, n = f"{name}_{n}", n + 1
            row = {"name": unique, "path": path, "kind": "experimental", "reference_standards": False,
                   "rt_unit": "min", "rt_tolerance_s": 10.0}
            table = pd.concat([table, pd.DataFrame([row])], ignore_index=True)
        self.libraries.value = table
        self._say("Library added (not saved yet): check its kind and whether it holds reference standards.", "info")

    def _model_dropped(self, event) -> None:
        self.config.scoring.ms2deepscore_model = event.new[-1]
        self._say(f"MS2DeepScore model set to {event.new[-1]} (not saved yet).", "info")

    def _ms2query_dropped(self, event) -> None:
        self.config.ms2query.models_dir = str(MS2QUERY_DROP_DIR)
        self._say(f"MS2Query will use the files in {MS2QUERY_DROP_DIR} (not saved yet).", "info")

    def _files_rejected(self, event) -> None:
        self._say(f"Not added (wrong file type): {', '.join(event.new[len(event.old):])}", "warning")

    def _read_library_summary(self) -> pd.DataFrame:
        path = self.project.results_dir / "annotations" / "library_summary.tsv"
        return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame({"library": []})

    # ---- SIRIUS account ----

    def _sirius_status_text(self) -> str:
        account = load_sirius_account(self.account_file)
        return f"Saved: {account['username']}." if account else "No account saved."

    def save_sirius_login(self) -> bool:
        if not (self.sirius_user.value.strip() and self.sirius_password.value):
            self.sirius_status.object = "Enter the account e-mail and password first."
            return False
        if not self.sirius_terms.value:
            self.sirius_status.object = "SIRIUS can only log in once its terms of service are accepted."
            return False
        save_sirius_account(self.sirius_user.value.strip(), self.sirius_password.value, True, self.account_file)
        self.sirius_password.value = ""
        self.sirius_password.placeholder = "saved"
        self.sirius_status.object = self._sirius_status_text()
        return True

    def forget_sirius_login(self) -> None:
        forget_sirius_account(self.account_file)
        self.sirius_password.placeholder = ""
        self.sirius_status.object = self._sirius_status_text()

    # ---- libraries ----

    def _add_library(self, event) -> None:
        row = pd.DataFrame([{"name": f"library_{len(self.libraries.value) + 1}", "path": "", "kind": "experimental",
                             "reference_standards": False, "rt_unit": "min", "rt_tolerance_s": 10.0}])
        self.libraries.value = pd.concat([self.libraries.value, row], ignore_index=True)

    def _remove_libraries(self, event) -> None:
        self.libraries.value = self.libraries.value.drop(index=self.libraries.selection).reset_index(drop=True)
        self.libraries.selection = []

    def _library_entries(self) -> list[dict]:
        """The library table as settings entries; the files must exist."""
        entries = []
        for row in self.libraries.value.to_dict("records"):
            path = Path(str(row["path"])).expanduser()
            if not path.is_file():
                raise ValueError(f"Library '{row['name']}': file not found: '{row['path']}'")
            entries.append({"name": str(row["name"]).strip(), "path": str(path.resolve()), "kind": row["kind"],
                            "reference_standards": bool(row["reference_standards"]), "rt_unit": row["rt_unit"],
                            "rt_tolerance_s": float(row["rt_tolerance_s"])})
        return entries

    def save(self) -> bool:
        """Validate and write samples.tsv and project.yaml. Returns False on errors."""
        try:
            self.config.adducts = AdductSettings(
                polarity=self.polarity.value,
                adducts=[{"name": r["name"], "openms": r["openms"], "probability": float(r["probability"])}
                         for r in self.adducts.value.to_dict("records")],
            )
            # The section object stays (its parameter card is bound to it): set and check the list.
            self.config.library_search.libraries = self._library_entries()
            self.config.library_search.check_libraries()
            self.project.save_samples(self.samples.value)
            self.project.save_config(self.config)
        except (ValueError, KeyError) as error:
            self._say(f"Not saved: {error}", "danger")
            return False
        self._say("Samples and parameters saved.")
        return True
