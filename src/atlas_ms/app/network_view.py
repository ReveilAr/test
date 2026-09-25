"""
The Network tab: molecular network, feature list and the plots of the
selected feature, all linked through one selection.

* Clicking (or box-selecting) nodes selects features; so does clicking rows
  of the feature list. Both views show the same selection.
* The plots show the first selected feature: its MS2 spectrum (the one used
  for networking) and its extracted ion chromatogram in every run.

Node colours are computed here, as one colour per node, which keeps the
plotting code free of colour-mapping machinery.
"""

import holoviews as hv
import numpy as np
import pandas as pd
import panel as pn
import param
from bokeh.models import Span
from bokeh.palettes import Category20, Viridis256

from atlas_ms.app.data import ChromatogramReader, Results, log_intensity

hv.extension("bokeh")

GREY = "#c8c8c8"
EDGE_COLORS = {"spectral": "#8c8c8c", "ion_identity": "#e6550d"}
COLOR_OPTIONS = ["family", "community", "intensity", "gap-filled values", "spectrum QC"]
TABLE_COLUMNS = ["feature_id", "mz", "rt", "intensity", "ion", "family", "community", "n_ms2", "n_gap_filled", "qc"]


def group_colors(groups: pd.Series) -> list[str]:
    """One colour per group number (cycling through 20 colours); grey for -1 (no group)."""
    palette = Category20[20]
    return [GREY if g < 1 else palette[(int(g) - 1) % len(palette)] for g in groups]


def node_colors(nodes: pd.DataFrame, color_by: str) -> list[str]:
    """The colour of every node for the chosen "colour by" option."""
    if color_by in ("family", "community"):
        return group_colors(nodes[color_by])
    if color_by == "intensity":
        values = log_intensity(nodes["intensity"])
        low, high = values.min(), values.max()
        scaled = (values - low) / (high - low) if high > low else values * 0
        return [Viridis256[int(v * 255)] for v in scaled]
    if color_by == "gap-filled values":
        return ["#e6550d" if n > 0 else "#3182bd" for n in nodes["n_gap_filled"]]
    # spectrum QC: red when the spectrum was kept out of spectral edges
    return ["#3182bd" if qc == "" else "#de2d26" for qc in nodes["qc"].fillna("")]


class NetworkTab(param.Parameterized):
    """State and views of the Network tab."""

    color_by = param.Selector(default="family", objects=COLOR_OPTIONS, doc="Node colour")
    selected = param.List(default=[], doc="Selected feature ids")

    def __init__(self, results: Results, reader: ChromatogramReader, **params):
        super().__init__(**params)
        self.results = results
        self.reader = reader
        # Network nodes: the features with MS2 (they have layout coordinates).
        self.nodes = results.features.dropna(subset=["x"]).reset_index(drop=True)
        # Node size grows with log intensity, from 8 to 20 pixels.
        logs = log_intensity(self.nodes["intensity"])
        self.nodes["size"] = np.interp(logs, (logs.min(), max(logs.max(), logs.min() + 1e-9)), (8, 20))
        self.table = pn.widgets.Tabulator(
            table_view(results.features),
            selectable=True, disabled=True, show_index=False, pagination="local", page_size=12,
            sizing_mode="stretch_width", height=420,
        )
        self.table.param.watch(self._table_selection, "selection")

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
            "score": edges["score"].to_numpy(),
        })
        return hv.Segments(data, ["x0", "y0", "x1", "y1"], ["color", "score"]).opts(
            color="color", line_width=1.5, alpha=0.8)

    def _points(self, color_by: str) -> hv.Points:
        data = self.nodes.assign(color=node_colors(self.nodes, color_by), ion=self.nodes["ion"].fillna(""))
        return hv.Points(data, ["x", "y"], ["feature_id", "mz", "rt", "family", "ion", "color", "size"]).opts(
            color="color", size="size", line_color="#444444", line_width=0.5,
            tools=["tap", "box_select", "hover"], nonselection_alpha=0.6,
        )

    def _highlight(self, selected: list) -> hv.Points:
        """Rings around the selected nodes (whichever view selected them)."""
        chosen = self.nodes[self.nodes["feature_id"].isin(selected)]
        return hv.Points(chosen, ["x", "y"]).opts(
            size=22, fill_alpha=0, line_color="#de2d26", line_width=3)

    def network(self) -> pn.pane.HoloViews:
        points = hv.DynamicMap(pn.bind(self._points, self.param.color_by))
        # Clicks and box selections on the nodes -> selected feature ids.
        selection = hv.streams.Selection1D(source=points)
        selection.param.watch(self._network_selection, "index")
        highlight = hv.DynamicMap(pn.bind(self._highlight, self.param.selected))
        overlay = (self._edges() * points * highlight).opts(
            hv.opts.Overlay(xaxis=None, yaxis=None, responsive=True, height=620,
                            active_tools=["wheel_zoom"], title="Molecular network"))
        return pn.pane.HoloViews(overlay, sizing_mode="stretch_width")

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

    # ---- plots of the selected feature ----

    @param.depends("selected")
    def details(self):
        if not self.selected:
            return pn.pane.Markdown("Select a node or a row to see its spectrum and chromatograms.")
        feature_id = self.selected[0]
        feature = self.results.features.set_index("feature_id").loc[feature_id]
        text = f"**Feature {feature_id}**: m/z {feature['mz']:.4f}, RT {feature['rt']:.1f} s"
        if isinstance(feature.get("ion"), str):
            text += f", {feature['ion']}"
        if pd.notna(feature.get("family")) and feature["family"] > 0:
            text += f", family {int(feature['family'])}"
        plots = [pn.pane.Markdown(text)]
        if feature_id in self.results.spectra:
            mz, intensity = self.results.spectra[feature_id]
            spectrum = pd.DataFrame({"mz": mz, "intensity": intensity})
            figure = spectrum.plot(kind="spectrum", x="mz", y="intensity", backend="ms_bokeh",
                                   show_plot=False, title="MS2 spectrum")
            plots.append(_fit(figure))
        try:
            xic = self.reader.xic(feature["mz"], feature["rt"])
        except (OSError, RuntimeError) as error:  # e.g. raw files moved since processing
            plots.append(pn.pane.Alert(f"Chromatograms unavailable: {error}", alert_type="warning"))
            return pn.Column(*plots, sizing_mode="stretch_width")
        figure = xic.plot(kind="chromatogram", x="rt", y="intensity", by="sample", backend="ms_bokeh",
                          show_plot=False, title="Chromatograms (aligned RT)")
        figure.add_layout(Span(location=feature["rt"], dimension="height", line_dash="dashed", line_color="#888888"))
        plots.append(_fit(figure))
        return pn.Column(*plots, sizing_mode="stretch_width")

    def view(self) -> pn.viewable.Viewable:
        controls = pn.Row(pn.widgets.Select.from_param(self.param.color_by, width=200))
        return pn.Column(
            controls,
            pn.Row(self.network(), sizing_mode="stretch_width"),
            pn.Row(pn.Column(self.table, sizing_mode="stretch_width"),
                   pn.Column(self.details, sizing_mode="stretch_width"),
                   sizing_mode="stretch_width"),
            sizing_mode="stretch_width",
        )


def table_view(features: pd.DataFrame) -> pd.DataFrame:
    """The feature list as displayed: rounded values, whole numbers without decimals."""
    table = features[[c for c in TABLE_COLUMNS if c in features]].copy()
    table["mz"] = table["mz"].round(4)
    table["rt"] = table["rt"].round(1)
    table["intensity"] = table["intensity"].round(0)
    for column in ("family", "community"):
        if column in table:
            # Features without MS2 are not in the network: no family ("Int64" allows missing).
            table[column] = table[column].astype("Int64")
    return table


def _fit(figure) -> pn.pane.Bokeh:
    """A Bokeh figure stretched to the width of its container."""
    figure.sizing_mode = "stretch_width"
    figure.height = 280
    return pn.pane.Bokeh(figure, sizing_mode="stretch_width")
