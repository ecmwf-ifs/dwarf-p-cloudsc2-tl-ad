# -*- coding: utf-8 -*-
from __future__ import annotations
import dataclasses
from functools import partial
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from typing import TYPE_CHECKING

from scripts import plot_utils

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Optional


@dataclasses.dataclass
class Line:
    column: str
    color: str
    linestyle: str
    marker: str
    linewidth: float
    markersize: float
    markerfacecolor: Optional[str] = None
    markeredgecolor: Optional[str] = None
    label: Optional[str] = None

    def __post_init__(self):
        self.markerfacecolor = self.markerfacecolor or self.color
        self.markeredgecolor = self.markeredgecolor or self.color
        self.label = self.label or self.column


DEFAULT_LINEWIDTH = 1.5
DEFAULT_MARKERSIZE = 10
DefaultLine = partial(Line, linewidth=DEFAULT_LINEWIDTH, markersize=DEFAULT_MARKERSIZE, label=None)


@dataclasses.dataclass
class LinePool:
    csv_file_generator: Callable[[int], str]
    title: str
    lines: list[Line]
    reference_column: str = "fortran"


def get_csv_file(nx: int, mode: str) -> str:
    return f"20220805_dom_{mode}_{nx}_137_24.csv"


lines = [
    DefaultLine("gt:cpu_kfirst", "cyan", "-", ">"),
    DefaultLine("gt:gpu", "green", "-", "<"),
    DefaultLine("cuda", "purple", "-", "o"),
    DefaultLine("dace:gpu", "orange", "-", "s"),
]
line_pool_nl = LinePool(partial(get_csv_file, mode="nl"), "$\\mathbf{(a)}$ Non-linear", lines)
line_pool_tl = LinePool(
    partial(get_csv_file, mode="tl"), "$\\mathbf{(b)}$ Tangent-linear", lines[0:1] + lines[2:]
)
line_pool_ad = LinePool(
    partial(get_csv_file, mode="ad"), "$\\mathbf{(c)}$ Adjoint", lines[0:1] + lines[2:]
)
nx_l = tuple(2**i for i in range(10, 18))


figure_properties = {
    "figsize": [14, 6],
    "fontsize": 16,
    "tight_layout": True,
}
axes_properties = {
    "fontsize": 16,
    "x_label": "Number of columns [-]",
    "x_scale": "log",
    "x_lim": None,
    "x_ticks": nx_l,
    "x_tick_labels": tuple(f"$2^{{{int(np.log2(nx))}}}$" for nx in nx_l),
    "y_label": "Speed-up w.r.t. Fortran [-]",
    "y_lim": [0, 18],
    "y_ticks": [0, 1, 4, 8, 12, 16, 20],
    "y_tick_labels": [0, 1, 4, 8, 12, 16, 20],
    "legend_on": True,
    "legend_fontsize": 14,
    "legend_framealpha": 1.0,
    "legend_loc": "upper left",
    "legend_ncol": 1,
    "grid_on": True,
    "grid_properties": {"linestyle": ":"},
}


def plot_lines(ax: plt.Axes, line_pool: LinePool) -> dict[str, list[float]]:
    speedups = {line.column: [] for line in line_pool.lines}

    # collect data
    for nx in nx_l:
        csv_file = line_pool.csv_file_generator(nx)
        df = pd.read_csv(os.path.join("timings", csv_file))
        data = df[df[line_pool.reference_column].notna()][line_pool.reference_column]
        reference_runtime = sum(data) / len(data)

        for line in line_pool.lines:
            if line.column in df.columns:
                data = df[df[line.column].notna()][line.column]
                runtime = sum(data) / len(data)
                speedups[line.column].append(reference_runtime / runtime)
            else:
                speedups[line.column].append(math.nan)

    # plot
    ax.plot(nx_l, (1,) * len(nx_l), "k-", linewidth=1)
    for line in line_pool.lines:
        ax.plot(
            nx_l,
            speedups[line.column],
            line.color,
            linestyle=line.linestyle,
            linewidth=line.linewidth,
            marker=line.marker,
            markersize=line.markersize,
            markerfacecolor=line.markerfacecolor,
            markeredgecolor=line.markeredgecolor,
            label=line.label,
        )

    return speedups


def fill_ax(ax: plt.Axes, line_pool: LinePool, axes_properties: dict) -> None:
    plot_lines(ax, line_pool)
    axes_properties["title_center"] = line_pool.title
    plot_utils.set_axes_properties(ax, **axes_properties)


def main():
    fig, ax = plot_utils.get_figure_and_axes(nrows=1, ncols=3, index=1, **figure_properties)
    fill_ax(ax, line_pool_nl, axes_properties)

    axes_properties["y_label"] = ""
    axes_properties["y_tick_labels"] = ()
    axes_properties["legend_on"] = False

    _, ax = plot_utils.get_figure_and_axes(fig=fig, nrows=1, ncols=3, index=2, **figure_properties)
    fill_ax(ax, line_pool_tl, axes_properties)

    _, ax = plot_utils.get_figure_and_axes(fig=fig, nrows=1, ncols=3, index=3, **figure_properties)
    fill_ax(ax, line_pool_ad, axes_properties)

    plot_utils.set_figure_properties(fig, **figure_properties)

    plt.show()


if __name__ == "__main__":
    main()
