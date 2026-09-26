"""
The Annotation tab: a closer look at the feature selected in the Network
tab (the selection is shared).

* Header: the feature, its best annotation (Schymanski level, and the
  structural level for lipids), its flags, and its family's class consensus.
* Candidates: every candidate of every source, most confident first, with
  its final level, the level its source proposed, scores and evidence.
* Evidence plot for the candidate clicked in the list:
  - library hit: mirror plot, the feature's spectrum against the library's;
  - lipid rule: the feature's spectrum with the diagnostic ions that were
    found marked (fragments at their m/z, neutral losses at precursor - loss).
* Family: the other members of the feature's molecular family and their
  annotations; clicking one selects it everywhere.
"""

import json

import numpy as np
import pandas as pd
import panel as pn

from atlas_ms.annotation.schema import LEVEL_RANK
from atlas_ms.app.data import Results
from atlas_ms.app.plots import mirror_plot, spectrum_plot

CANDIDATE_COLUMNS = ["level", "name", "source", "score", "matched_peaks", "mz_error_ppm", "rt_error_s",
                     "adduct", "formula", "lipid_class", "lipid_level", "proposed_level", "evidence_text"]
CANDIDATE_TITLES = {"matched_peaks": "matched ions", "mz_error_ppm": "m/z error (ppm)", "rt_error_s": "RT error (s)",
                    "lipid_class": "lipid class", "lipid_level": "structural level",
                    "proposed_level": "proposed", "evidence_text": "evidence"}
# Numbers shown with a fixed number of decimals; blank when not applicable.
DECIMALS = {"score": 2, "matched_peaks": 0, "mz_error_ppm": 1, "rt_error_s": 0}
FAMILY_COLUMNS = ["feature_id", "mz", "rt_min", "annotation", "level", "lipid_class", "annotation_flags"]
FAMILY_TITLES = {"feature_id": "feature", "mz": "m/z", "rt_min": "RT (min)", "lipid_class": "lipid class",
                 "annotation_flags": "flags"}


def evidence_text(evidence: str) -> str:
    """The JSON evidence of a candidate as one readable line."""
    parts = []
    for key, value in (json.loads(evidence) if evidence else {}).items():
        if isinstance(value, list):
            value = ", ".join(str(v) for v in value) or "-"
        parts.append(f"{key.replace('_', ' ')}: {value}")
    return "; ".join(parts)


def evidence_marks(evidence: str, precursor: float) -> list[tuple[float, str]]:
    """
    The ions a lipid rule found, as (m/z, text) marks: fragments at their
    m/z ("184.0733"), neutral losses at precursor - loss ("NL 141.0191").
    Chain evidence is text only (it names several ions).
    """
    marks = []
    for item in (json.loads(evidence) if evidence else {}).get("found", []):
        if item.startswith("NL "):
            marks.append((precursor - float(item[3:]), item))
        else:
            try:
                marks.append((float(item), item))
            except ValueError:
                pass
    return marks


def candidate_view(rows: pd.DataFrame) -> pd.DataFrame:
    """The candidate list as displayed: rounded numbers, blanks instead of NaN."""
    table = rows[CANDIDATE_COLUMNS].copy()
    for column, decimals in DECIMALS.items():
        table[column] = [("" if pd.isna(v) else f"{v:.{decimals}f}") for v in table[column]]
    return table


def has_peaks(values) -> bool:
    """A reference spectrum column holds a list (or array) of m/z, or nothing."""
    return isinstance(values, (list, np.ndarray)) and len(values) > 0


class AnnotationTab:
    """State and views of the Annotation tab. ``network`` holds the shared selection."""

    def __init__(self, results: Results, network):
        self.results = results
        self.network = network
        self.features = results.features.set_index("feature_id")
        candidates = results.candidates
        if not candidates.empty:
            candidates = candidates.assign(
                evidence_text=candidates["evidence"].map(evidence_text),
                _rank=candidates["level"].map(LEVEL_RANK),
            ).sort_values(["feature_id", "_rank", "score"], ascending=[True, True, False])
        self.candidates = candidates
        self._rows = pd.DataFrame()  # candidates of the selected feature (all columns)

        self.table = pn.widgets.Tabulator(
            pd.DataFrame(columns=CANDIDATE_COLUMNS), selectable=1, disabled=True, show_index=False,
            titles=CANDIDATE_TITLES, sizing_mode="stretch_width", pagination="local", page_size=15,
        )
        self.family_table = pn.widgets.Tabulator(
            pd.DataFrame(columns=FAMILY_COLUMNS), selectable=1, disabled=True, show_index=False,
            titles=FAMILY_TITLES, sizing_mode="stretch_width", pagination="local", page_size=10,
        )
        self.family_table.param.watch(self._select_member, "selection")
        network.param.watch(self._show_feature, "selected")
        self._show_feature()

    # ---- selection ----

    def _feature_id(self) -> int | None:
        return self.network.selected[0] if self.network.selected else None

    def _show_feature(self, event=None) -> None:
        """Fill the candidate and family lists for the selected feature."""
        feature_id = self._feature_id()
        rows = pd.DataFrame()
        if feature_id is not None and not self.candidates.empty:
            rows = self.candidates[self.candidates["feature_id"] == feature_id].reset_index(drop=True)
        self._rows = rows
        self.table.value = candidate_view(rows) if not rows.empty else pd.DataFrame(columns=CANDIDATE_COLUMNS)
        self.table.selection = [0] if not rows.empty else []

        members = pd.DataFrame(columns=FAMILY_COLUMNS)
        if feature_id is not None and self.features.loc[feature_id, "family"] > 0:
            family = self.features[self.features["family"] == self.features.loc[feature_id, "family"]]
            members = family.reset_index().assign(rt_min=family["rt"].to_numpy() / 60.0)[FAMILY_COLUMNS]
            members = members.round({"mz": 4, "rt_min": 2}).sort_values("feature_id", ignore_index=True)
        self.family_table.value = members

    def _select_member(self, event) -> None:
        """Clicking a family member selects it (in every tab)."""
        if event.new:
            feature_id = int(self.family_table.value.iloc[event.new[0]]["feature_id"])
            if [feature_id] != self.network.selected:
                self.network.selected = [feature_id]

    # ---- views ----

    def header(self, selected) -> pn.pane.Markdown:
        feature_id = self._feature_id()
        if feature_id is None:
            return pn.pane.Markdown("Select a feature in the Network tab (node or table row).")
        feature = self.features.loc[feature_id]
        text = f"### Feature {feature_id}: m/z {feature['mz']:.4f}, RT {feature['rt'] / 60:.2f} min\n"
        text += f"**Best annotation:** {feature['annotation_label'] or 'none'}"
        if feature["annotation_flags"]:
            text += f"  \n**Flags:** {feature['annotation_flags'].replace(';', ', ')}"
        if feature["family_class"]:
            text += (f"  \n**Family {int(feature['family'])} class consensus:** {feature['family_class']} "
                     f"({feature['family_class_score']:.0%} of its annotated members)")
        return pn.pane.Markdown(text)

    def evidence(self, value, selection) -> pn.viewable.Viewable:
        """Plot of the evidence behind the clicked candidate."""
        if self._rows.empty or not selection:
            return pn.pane.Markdown("No annotation candidate for this feature." if self._feature_id() else "")
        candidate = self._rows.iloc[selection[0]]
        feature_id = int(candidate["feature_id"])
        spectrum = self.results.spectra.get(feature_id)
        if spectrum is None or len(spectrum[0]) == 0:
            return pn.pane.Markdown("This feature has no MS2 spectrum: the candidate matches its m/z only.")
        if has_peaks(candidate["reference_mz"]):
            reference = (np.asarray(candidate["reference_mz"], dtype=float),
                         np.asarray(candidate["reference_intensity"], dtype=float))
            figure = mirror_plot(spectrum, reference, title=f"Feature {feature_id} (top) vs {candidate['name']} "
                                                            f"({candidate['source']}, bottom)")
        else:
            marks = evidence_marks(candidate["evidence"], float(self.features.loc[feature_id, "mz"]))
            figure = spectrum_plot(*spectrum, title=f"Feature {feature_id}: evidence for {candidate['name']}",
                                   marks=marks)
        return pn.pane.Bokeh(figure, sizing_mode="stretch_width")

    def view(self) -> pn.viewable.Viewable:
        if self.candidates.empty:
            return pn.pane.Markdown("No annotations yet: they are made by the pipeline (press **Run**).")
        return pn.Column(
            pn.bind(self.header, self.network.param.selected),
            pn.pane.Markdown("#### Candidates (click one to see its evidence)"),
            self.table,
            pn.bind(self.evidence, self.table.param.value, self.table.param.selection),
            pn.pane.Markdown("#### Molecular family"),
            self.family_table,
            sizing_mode="stretch_width",
        )
