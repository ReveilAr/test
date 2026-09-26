"""
The Network tab: molecular network, plots of the selected features, and the
feature list, all linked through one selection.

* Selecting: click a node (shift-click to add more), box-select nodes, or
  click rows of the feature list (ctrl-click to add more). All views show
  the same selection, circled in red in the network.
* Plots (right of the network): the MS2 spectrum of the first selected
  feature, or a mirror plot of the first two when two are selected; and the
  first feature's chromatogram in every run.
* Fragment search: features whose MS2 contains given fragments (or neutral
  losses) are circled in green.
* The feature list has a search box above every column.

Node colours and sizes are computed here, one value per node, which keeps
the plotting code free of colour-mapping machinery. Retention times are
shown in minutes (the pipeline itself works in seconds).

Overlapping nodes: node sizes are radii in layout units (like the x and y
coordinates), not pixels. The layout keeps nodes at least MIN_DISTANCE
apart and no radius exceeds MIN_DISTANCE / 2, so two nodes never overlap,
at any zoom: the mouse is only ever over one node, and its tooltip is the
only one shown. Zooming in makes the nodes bigger on screen, as in Cytoscape.
"""

import re

import holoviews as hv
import numpy as np
import pandas as pd
import panel as pn
import param
from bokeh.models.widgets.tables import NumberFormatter
from bokeh.palettes import Category20, Viridis256

from atlas_ms.app.data import ChromatogramReader, Results, find_fragments, log_intensity
from atlas_ms.app.plots import chromatogram_plot, mirror_plot, spectrum_plot
from atlas_ms.network.graph import MIN_DISTANCE

hv.extension("bokeh")

GREY = "#c8c8c8"
EDGE_COLORS = {"spectral": "#8c8c8c", "ion_identity": "#e6550d"}
SELECTED_RING = "#de2d26"  # red
MATCH_RING = "#1b9e77"  # green
COLOR_OPTIONS = ["family", "community", "retention time", "intensity", "gap-filled values", "spectrum QC"]
SAME_SIZE, MEAN_SIZE = "same size", "mean intensity"
TABLE_COLUMNS = ["feature_id", "mz", "rt_min", "intensity", "ion", "family", "community",
                 "n_ms2", "n_gap_filled", "qc"]
TABLE_TITLES = {"feature_id": "feature", "mz": "m/z", "rt_min": "RT (min)", "n_ms2": "MS2 spectra",
                "n_gap_filled": "gap-filled", "qc": "spectrum QC"}
TABLE_FORMATS = {  # fixed decimals, no thousands separators: easier to read and to search
    "mz": NumberFormatter(format="0.0000", text_align="right"),
    "rt_min": NumberFormatter(format="0.00", text_align="right"),
    "intensity": NumberFormatter(format="0", text_align="right"),
}
HOVER = [("feature", "@feature_id"), ("m/z", "@mz{0.0000}"), ("RT (min)", "@rt_min{0.00}"),
         ("family", "@family"), ("ion", "@ion")]


# --------------------------------------------------------------------------
# Colours and sizes
# --------------------------------------------------------------------------

def group_colors(groups: pd.Series) -> list[str]:
    """One colour per group number (cycling through 20 colours); grey for -1 (no group)."""
    palette = Category20[20]
    return [GREY if g < 1 else palette[(int(g) - 1) % len(palette)] for g in groups]


def gradient_colors(values: pd.Series) -> list[str]:
    """Viridis colours from the lowest (dark blue) to the highest value (yellow)."""
    low, high = values.min(), values.max()
    scaled = (values - low) / (high - low) if high > low else values * 0
    return [Viridis256[int(v * 255)] for v in scaled]


def node_colors(nodes: pd.DataFrame, color_by: str) -> list[str]:
    """The colour of every node for the chosen "colour by" option."""
    if color_by in ("family", "community"):
        return group_colors(nodes[color_by])
    if color_by == "retention time":
        return gradient_colors(nodes["rt"])
    if color_by == "intensity":
        return gradient_colors(log_intensity(nodes["intensity"]))
    if color_by == "gap-filled values":
        return ["#e6550d" if n > 0 else "#3182bd" for n in nodes["n_gap_filled"]]
    # spectrum QC: red when the spectrum was kept out of spectral edges
    return ["#3182bd" if qc == "" else "#de2d26" for qc in nodes["qc"].fillna("")]


# Node radii in layout units. The largest is just under half the minimum
# distance between two nodes, so that nodes never overlap.
MAX_RADIUS = MIN_DISTANCE / 2 - 0.05
MIN_RADIUS = 0.45 * MAX_RADIUS
ABSENT_RADIUS = 0.25 * MAX_RADIUS  # feature not detected in the chosen sample


def node_radii(intensities: pd.Series) -> np.ndarray:
    """
    Node radii growing with log intensity, from MIN_RADIUS to MAX_RADIUS.
    Features absent from the chosen sample get the smallest radius.
    """
    logs = np.log10(intensities.where(intensities > 0))
    low, high = np.nanmin(logs), np.nanmax(logs)
    radii = np.interp(logs, (low, max(high, low + 1e-9)), (MIN_RADIUS, MAX_RADIUS))
    return np.where(np.isnan(logs), ABSENT_RADIUS, radii)


def table_view(features: pd.DataFrame) -> pd.DataFrame:
    """The feature list as displayed: RT in minutes, rounded values, whole numbers without decimals."""
    table = features.assign(rt_min=features["rt"] / 60.0)
    table = table[[c for c in TABLE_COLUMNS if c in table]].copy()
    table["mz"] = table["mz"].round(4)
    table["rt_min"] = table["rt_min"].round(2)
    table["intensity"] = table["intensity"].round(0)
    for column in ("family", "community"):
        if column in table:
            # Features without MS2 are not in the network: no family ("Int64" allows missing).
            table[column] = table[column].astype("Int64")
    return table


def parse_masses(text: str) -> list[float]:
    """m/z values typed as "184.0733, 104.107" (commas, semicolons or spaces)."""
    return [float(value) for value in re.split(r"[,;\s]+", text.strip()) if value]


# --------------------------------------------------------------------------
# The tab
# --------------------------------------------------------------------------

class NetworkTab(param.Parameterized):
    """State and views of the Network tab."""

    color_by = param.Selector(default="family", objects=COLOR_OPTIONS, doc="Node colour")
    size_by = param.Selector(default=MEAN_SIZE, objects=[SAME_SIZE, MEAN_SIZE], doc="Node size")
    mirror = param.Boolean(default=True, doc="Mirror plot when two features are selected")
    selected = param.List(default=[], doc="Selected feature ids")
    matched = param.List(default=[], doc="Feature ids found by the fragment search")

    def __init__(self, results: Results, reader: ChromatogramReader, fragment_tolerance: float = 0.01, **params):
        super().__init__(**params)
        self.results = results
        self.reader = reader
        # Node size can follow the feature's area in any single sample.
        self.param.size_by.objects = [SAME_SIZE, MEAN_SIZE] + list(results.quant.columns)

        # Network nodes: the features with MS2 (they have layout coordinates).
        nodes = results.features.dropna(subset=["x"]).reset_index(drop=True)
        self.nodes = nodes.assign(
            rt_min=nodes["rt"] / 60.0,
            ion=nodes["ion"].fillna(""),
            family=nodes["family"].astype(int),
            community=nodes["community"].astype(int),
        )

        # ---- fragment search ----
        self.fragment_text = pn.widgets.TextInput(label="MS2 search (m/z, comma-separated)",
                                                  placeholder="184.0733, 104.1070", width=260)
        self.fragment_kind = pn.widgets.RadioButtonGroup(options=["fragments", "neutral losses"], value="fragments")
        self.fragment_tolerance = pn.widgets.FloatInput(label="± Da", value=fragment_tolerance, step=0.005,
                                                        start=0.0001, width=90)
        search = pn.widgets.Button(label="Search", color="primary", width=90)
        clear = pn.widgets.Button(label="Clear", width=80)
        search.on_click(lambda event: self.search_fragments())
        clear.on_click(lambda event: self.clear_search())
        self.search_message = pn.pane.Markdown("")

        # ---- feature list ----
        table = table_view(results.features)
        self.table = pn.widgets.Tabulator(
            table, selectable=True, disabled=True, show_index=False, titles=TABLE_TITLES, formatters=TABLE_FORMATS,
            pagination="local", page_size=25, sizing_mode="stretch_width",
            # A search box above every column; "like" = contains, also for numbers (e.g. "760.58").
            header_filters={c: {"type": "input", "func": "like", "placeholder": "search"} for c in table.columns},
        )
        self.table.param.watch(self._table_selection, "selection")

        self.controls = pn.Row(
            pn.widgets.Select.from_param(self.param.color_by, label="Colour by", width=170),
            pn.widgets.Select.from_param(self.param.size_by, label="Size by", width=170),
            self.fragment_text, self.fragment_tolerance,
            pn.Column(self.fragment_kind, pn.Row(search, clear)),
            sizing_mode="stretch_width",
        )

    # ---- network ----

    def _edges(self) -> hv.Segments:
        position = self.nodes.set_index("feature_id")[["x", "y"]]
        edges = self.results.edges[self.results.edges["source"].isin(position.index)
                                   & self.results.edges["target"].isin(position.index)]
        data = pd.DataFrame({
            "x0": position.loc[edges["source"], "x"].to_numpy(),
            "y0": position.loc[edges["source"], "y"].to_numpy(),
            "x1": position.loc[edges["target"], "x"].to_numpy(),
            "y1": position.loc[edges["target"], "y"].to_numpy(),
            "color": [EDGE_COLORS[t] for t in edges["edge_type"]],
        })
        return hv.Segments(data, ["x0", "y0", "x1", "y1"], ["color"]).opts(
            color="color", line_width=1, alpha=0.8)

    def _radii(self, size_by: str) -> np.ndarray:
        if size_by == SAME_SIZE:
            return np.full(len(self.nodes), 0.8 * MAX_RADIUS)
        if size_by == MEAN_SIZE:
            return node_radii(self.nodes["intensity"])
        return node_radii(self.results.quant[size_by].reindex(self.nodes["feature_id"]).reset_index(drop=True))

    def _points(self, color_by: str, size_by: str) -> hv.Points:
        data = self.nodes.assign(color=node_colors(self.nodes, color_by), radius=self._radii(size_by))
        return hv.Points(data, ["x", "y"], ["feature_id", "mz", "rt_min", "family", "ion", "color", "radius"]).opts(
            color="color", radius="radius", line_color="#444444", line_width=0.5,
            tools=["tap", "box_select", "hover"], hover_tooltips=HOVER, nonselection_alpha=0.6,
        )

    def _rings(self, feature_ids: list, color: str, radius: float) -> hv.Points:
        """
        Rings around some nodes (the selection, or the fragment-search
        matches). Their thick line keeps them visible when zoomed out.
        """
        chosen = self.nodes[self.nodes["feature_id"].isin(feature_ids)]
        return hv.Points(chosen, ["x", "y"]).opts(radius=radius, fill_alpha=0, line_color=color, line_width=3)

    def network(self) -> pn.pane.HoloViews:
        points = hv.DynamicMap(pn.bind(self._points, self.param.color_by, self.param.size_by))
        # Clicks and box selections on the nodes -> selected feature ids.
        selection = hv.streams.Selection1D(source=points)
        selection.param.watch(self._network_selection, "index")
        matches = hv.DynamicMap(pn.bind(self._rings, self.param.matched, MATCH_RING, 1.2 * MAX_RADIUS))
        selected = hv.DynamicMap(pn.bind(self._rings, self.param.selected, SELECTED_RING, 1.5 * MAX_RADIUS))
        overlay = (self._edges() * points * matches * selected).opts(
            # data_aspect=1 keeps the layout undistorted (round nodes, true distances).
            hv.opts.Overlay(xaxis=None, yaxis=None, responsive=True, min_height=650, data_aspect=1,
                            active_tools=["wheel_zoom"], title="Molecular network (scroll to zoom)"))
        return pn.pane.HoloViews(overlay, sizing_mode="stretch_width", min_height=650)

    # ---- fragment search ----

    def search_fragments(self) -> None:
        try:
            masses = parse_masses(self.fragment_text.value)
        except ValueError:
            self.search_message.object = "Type m/z values separated by commas, e.g. `184.0733, 104.107`."
            return
        if not masses:
            return self.clear_search()
        losses = self.fragment_kind.value == "neutral losses"
        precursors = dict(zip(self.nodes["feature_id"], self.nodes["mz"])) if losses else None
        spectra = {f: s for f, s in self.results.spectra.items() if precursors is None or f in precursors}
        self.matched = find_fragments(spectra, masses, self.fragment_tolerance.value, precursors)
        what = "neutral loss(es)" if losses else "fragment(s)"
        self.search_message.object = (f"**{len(self.matched)} features** have all {len(masses)} {what} "
                                      f"(circled in green).")

    def clear_search(self) -> None:
        self.matched = []
        self.fragment_text.value = ""
        self.search_message.object = ""

    # ---- selection, kept in sync between network and table ----

    def _network_selection(self, event) -> None:
        self.selected = self.nodes.loc[event.new, "feature_id"].astype(int).tolist()

    def _table_selection(self, event) -> None:
        ids = self.table.value.iloc[event.new]["feature_id"].astype(int).tolist()
        if ids != self.selected:
            self.selected = ids

    @param.depends("selected", watch=True)
    def _select_rows(self) -> None:
        rows = np.flatnonzero(self.table.value["feature_id"].isin(self.selected)).tolist()
        if rows != self.table.selection:
            self.table.selection = rows

    # ---- plots of the selected features ----

    def _edge_text(self, first: int, second: int) -> str:
        """Similarity of two features if the network links them."""
        edges = self.results.edges
        link = edges[((edges["source"] == first) & (edges["target"] == second))
                     | ((edges["source"] == second) & (edges["target"] == first))]
        if link.empty:
            return "not linked in the network"
        edge = link.iloc[0]
        if edge["edge_type"] == "ion_identity":
            return "linked as adducts of the same molecule"
        matched = f", {int(edge['matched_peaks'])} matched fragments" if pd.notna(edge["matched_peaks"]) else ""
        return f"similarity {edge['score']:.2f}{matched}"

    @param.depends("selected", "mirror")
    def details(self):
        if not self.selected:
            return pn.pane.Markdown("Click a node or a row to see its spectrum and chromatograms. "
                                    "Select two (shift-click / ctrl-click) for a mirror plot.")
        features = self.results.features.set_index("feature_id")
        feature_id = self.selected[0]
        feature = features.loc[feature_id]
        text = f"**Feature {feature_id}**: m/z {feature['mz']:.4f}, RT {feature['rt'] / 60:.2f} min"
        if isinstance(feature.get("ion"), str):
            text += f", {feature['ion']}"
        if pd.notna(feature.get("family")) and feature["family"] > 0:
            text += f", family {int(feature['family'])}"
        plots = [pn.pane.Markdown(text)]

        spectra = self.results.spectra
        if self.mirror and len(self.selected) >= 2 and all(f in spectra for f in self.selected[:2]):
            other = self.selected[1]
            plots.append(pn.pane.Markdown(
                f"Top: feature {feature_id}. Bottom: feature {other} (m/z {features.loc[other, 'mz']:.4f}): "
                f"{self._edge_text(feature_id, other)}."))
            figure = mirror_plot(spectra[feature_id], spectra[other], title="MS2 mirror plot")
            plots.append(pn.pane.Bokeh(figure, sizing_mode="stretch_width"))
        elif feature_id in spectra:
            plots.append(pn.pane.Bokeh(spectrum_plot(*spectra[feature_id]), sizing_mode="stretch_width"))

        try:
            xic = self.reader.xic(feature["mz"], feature["rt"])
        except (OSError, RuntimeError) as error:  # e.g. raw files moved since processing
            plots.append(pn.pane.Alert(f"Chromatograms unavailable: {error}", alert_type="warning"))
        else:
            if xic.empty:  # no MS1 scan within the RT window (feature at the very end of the runs)
                plots.append(pn.pane.Alert("No MS1 scans around this retention time.", alert_type="warning"))
            else:
                plots.append(pn.pane.Bokeh(chromatogram_plot(xic, feature["rt"]), sizing_mode="stretch_width"))
        return pn.Column(*plots, sizing_mode="stretch_width")

    def view(self) -> pn.viewable.Viewable:
        mirror = pn.widgets.Checkbox.from_param(self.param.mirror, label="Mirror plot for two selected features")
        return pn.Column(
            self.controls,
            self.search_message,
            pn.Row(
                pn.Column(self.network(), sizing_mode="stretch_width"),
                pn.Column(mirror, self.details, width=560),
                sizing_mode="stretch_width",
            ),
            pn.pane.Markdown("### Features"),
            self.table,
            sizing_mode="stretch_width",
        )
