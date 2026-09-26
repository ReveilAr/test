"""
The Setup tab: sample metadata and processing parameters of the open project.

* Sample table: edit sample types (sample / blank / standard) and free
  metadata columns (saved as ATTRIBUTE_<name>, the GNPS convention).
* Presets: copy an instrument or adduct preset into the parameters.
* Adducts: polarity and the editable adduct list.
* Spectral libraries: the MSP / MGF / JSON files to search, with their kind
  (experimental or in-silico) and whether they are reference standards
  measured on the same method (the only way to reach level 1).
* Parameters: every section of ``atlas_ms.config``, as input boxes.

Nothing is written until "Save" (or Run, which saves first). The files are
validated by the same code as the pipeline, so an invalid value is reported
here rather than failing a run.
"""

import re
from pathlib import Path

import pandas as pd
import panel as pn
import param

from atlas_ms.config import LIBRARY_KINDS, SECTIONS, AdductSettings, load_presets
from atlas_ms.project import SAMPLE_TYPES, Project

# Sections shown as plain parameter boxes (presets and adducts have their own widgets).
PARAMETER_SECTIONS = [name for name in SECTIONS if name not in ("presets", "adducts")]
# Card titles that differ from the capitalised section name.
TITLES = {"spectrum_qc": "Spectrum QC", "lipids": "Lipid rules", "sirius": "SIRIUS", "ms2query": "MS2Query"}
HIDDEN = {"name", "libraries"}  # parameters with their own editor (libraries: the library table)
LIBRARY_COLUMNS = ["name", "path", "kind", "reference_standards", "rt_unit", "rt_tolerance_s"]

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

    def __init__(self, project: Project):
        self.project = project
        self.config = project.load_config()
        self.message = pn.pane.Alert("", alert_type="light", visible=False, sizing_mode="stretch_width")

        # ---- sample table ----
        self.samples = pn.widgets.Tabulator(
            project.load_samples().reset_index(drop=True),
            show_index=False, sizing_mode="stretch_width", layout="fit_data_stretch",
            # Sample names and files are fixed; types come from a list; metadata is free text.
            editors={"sample": None, "file": None,
                     "sample_type": {"type": "list", "values": list(SAMPLE_TYPES)}},
        )
        self.new_column = pn.widgets.TextInput(placeholder="new metadata column, e.g. species", width=260)
        add_column = pn.widgets.Button(label="Add column", width=110)
        add_column.on_click(self._add_column)

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

        save = pn.widgets.Button(label="Save", color="primary", width=120)
        save.on_click(lambda event: self.save())

        # One card per parameter section, bound to the config objects: editing
        # a box changes the section directly, applying a preset updates the boxes.
        self._cards = pn.Column(*self._parameter_cards(), sizing_mode="stretch_width")
        self.layout = pn.Column(
            self.message,
            pn.pane.Markdown("### Samples"),
            self.samples,
            pn.Row(self.new_column, add_column),
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
            pn.pane.Markdown("### Parameters"),
            self._cards,
            save,
            sizing_mode="stretch_width",
        )

    def _parameter_cards(self) -> list:
        cards = []
        for name in PARAMETER_SECTIONS:
            section = getattr(self.config, name)
            parameters = [name for name in section.param if name not in HIDDEN]
            cards.append(pn.Card(
                pn.Param(section, parameters=parameters, widgets=widgets_for(section), show_name=False),
                title=TITLES.get(name, name.replace("_", " ").capitalize()), collapsed=True,
                sizing_mode="stretch_width",
            ))
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
