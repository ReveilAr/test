"""
Spectrum and chromatogram plots of the app (pyopenms-viz, Bokeh backend).

pyopenms-viz draws the plots. Two things are handled here on top of it:

* peak labels: pyopenms-viz labels the N most intense peaks, which overlap
  when intense peaks are close together. ``label_peaks`` chooses them
  instead: the most intense first, skipping any peak too close (in m/z)
  to a label already placed. The labels are small and drawn on both halves
  of a mirror plot;
* units: retention times in minutes, intensity ticks with one decimal
  (e.g. 6.0e+5).
"""

import math

import numpy as np
import pandas as pd
from bokeh.models import BasicTickFormatter, Label, Span

FONTS = {  # compact text for plots shown side by side
    "title_font_size": 12,
    "xaxis_label_font_size": 11,
    "yaxis_label_font_size": 11,
    "xaxis_tick_font_size": 10,
    "yaxis_tick_font_size": 10,
}
HEIGHT = 300
MARK_COLOR = "#c51b8a"  # the app's accent colour
# Minimum m/z distance between two labels, as a fraction of the plotted m/z
# range: a label ("760.5851" at 8 pt) is about 45 px wide, and a plot in the
# app's side column about 480 px, so labels closer than ~10% would overlap.
LABEL_GAP = 0.12


def label_peaks(mz: np.ndarray, intensity: np.ndarray, max_labels: int = 6,
                min_relative: float = 0.05, min_gap: float | None = None) -> list[int]:
    """
    Indices of the peaks to label: from the most intense down, at least
    ``min_relative`` of the base peak, and at least ``min_gap`` m/z away from
    every label already chosen (default: LABEL_GAP of the spectrum's m/z
    range, so labels never overlap). At most ``max_labels``.
    """
    if len(mz) == 0:
        return []
    if min_gap is None:
        min_gap = LABEL_GAP * max(np.ptp(mz), 50.0)
    chosen = []
    for index in np.argsort(intensity)[::-1]:
        if intensity[index] < min_relative * intensity.max() or len(chosen) == max_labels:
            break
        if all(abs(mz[index] - mz[other]) >= min_gap for other in chosen):
            chosen.append(int(index))
    return chosen


def _add_labels(figure, mz: np.ndarray, heights: np.ndarray, mz_range: float, below: bool = False) -> None:
    """
    m/z labels (4 decimals) above the chosen peaks, or below for the lower
    half of a mirror plot. ``mz_range`` is the plotted m/z range (both
    spectra of a mirror plot share the axis).
    """
    for index in label_peaks(mz, np.abs(heights), min_gap=LABEL_GAP * max(mz_range, 50.0)):
        figure.add_layout(Label(
            x=mz[index], y=heights[index], text=f"{mz[index]:.4f}", text_font_size="8pt",
            text_align="center", text_baseline="top" if below else "bottom", y_offset=-3 if below else 3,
        ))


def _finish(figure):
    figure.sizing_mode = "stretch_width"
    figure.height = HEIGHT
    return figure


def spectrum_plot(mz: np.ndarray, intensity: np.ndarray, title: str = "MS2 spectrum", marks=()):
    """
    One MS2 spectrum, absolute intensities. ``marks`` are (m/z, text) pairs
    drawn as dashed lines with vertical text, e.g. the diagnostic ions that
    support a lipid annotation.
    """
    peaks = pd.DataFrame({"mz": mz, "intensity": intensity})
    figure = peaks.plot(
        kind="spectrum", x="mz", y="intensity", backend="ms_bokeh", show_plot=False,
        annotate_top_n_peaks=0, title=title, xlabel="m/z", ylabel="Intensity", **FONTS,
    )
    figure.yaxis.formatter = BasicTickFormatter(precision=1)  # 6.0e+5
    _add_labels(figure, peaks["mz"].to_numpy(), peaks["intensity"].to_numpy(), np.ptp(mz) if len(mz) else 0.0)
    for position, text in marks:
        figure.add_layout(Span(location=position, dimension="height", line_dash="dashed",
                               line_color=MARK_COLOR, line_width=1.5))
        figure.add_layout(Label(x=position, y=0, angle=math.pi / 2, text=text, text_font_size="8pt",
                                text_color=MARK_COLOR, x_offset=-3, y_offset=4))
    return _finish(figure)


def mirror_plot(top: tuple[np.ndarray, np.ndarray], bottom: tuple[np.ndarray, np.ndarray], title: str):
    """
    Two MS2 spectra face to face, in relative intensity (% of each base peak),
    so that spectra of different intensities can be compared.
    """
    def relative(mz, intensity):
        return pd.DataFrame({"mz": mz, "intensity": 100.0 * intensity / intensity.max()})

    upper, lower = relative(*top), relative(*bottom)
    figure = upper.plot(
        kind="spectrum", x="mz", y="intensity", backend="ms_bokeh", show_plot=False,
        reference_spectrum=lower, mirror_spectrum=True, annotate_top_n_peaks=0,
        title=title, xlabel="m/z", ylabel="Relative intensity (%)", **FONTS,
    )
    figure.yaxis.formatter = BasicTickFormatter(precision=1, use_scientific=False)
    both = np.concatenate([upper["mz"], lower["mz"]])
    mz_range = np.ptp(both) if len(both) else 0.0
    _add_labels(figure, upper["mz"].to_numpy(), upper["intensity"].to_numpy(), mz_range)
    _add_labels(figure, lower["mz"].to_numpy(), -lower["intensity"].to_numpy(), mz_range, below=True)
    return _finish(figure)


def chromatogram_plot(xic: pd.DataFrame, feature_rt_s: float):
    """
    Extracted ion chromatograms of all runs (``msdata.RunReader.xic``), in
    minutes on the aligned axis, with a dashed line at the feature's RT.
    """
    data = xic.assign(rt=xic["rt"] / 60.0)
    figure = data.plot(
        kind="chromatogram", x="rt", y="intensity", by="sample", backend="ms_bokeh", show_plot=False,
        title="Chromatograms (aligned RT)", xlabel="Retention time (min)", ylabel="Intensity", **FONTS,
    )
    figure.yaxis.formatter = BasicTickFormatter(precision=1)
    figure.add_layout(Span(location=feature_rt_s / 60.0, dimension="height",
                           line_dash="dashed", line_color="#888888"))
    return _finish(figure)
