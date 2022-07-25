# -*- coding: utf-8 -*-
import dataclasses
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from typing import Optional, Sequence

from scripts import plot_utils


@dataclasses.dataclass
class Bar:
    csv_file: str
    column: str
    color: str
    label: str
    hatch: Optional[str] = None


@dataclasses.dataclass
class Pool:
    name: str
    bars: Sequence[Bar]


pool1 = Pool(
    "Non-linear",
    [
        # fortran vs gt4py
        # Bar(
        #     "dom_nl_24_20220220.csv",
        #     "fortran",
        #     "blue",
        #     "Fortran",
        # ),
        # Bar(
        #     "dom_nl_24_20220220.csv",
        #     "gtc:gt:cpu_kfirst",
        #     "cyan",
        #     "GT4Py CPU",
        # ),
        # Bar("dom_nl_24_20220220.csv", "gtc:gt:gpu", "green", "GT4Py GPU"),
        # single-core vs multi-core
        # Bar(
        #     "dom_nl_1_20220220.csv",
        #     "fortran",
        #     "blue",
        #     "Fortran (1 thread)",
        #     "//",
        # ),
        # Bar("dom_nl_24_20220220.csv", "fortran", "blue", "Fortran (24 threads)"),
        # Bar(
        #     "dom_nl_1_20220220.csv",
        #     "gtc:gt:cpu_kfirst",
        #     "cyan",
        #     "GT4Py CPU (1 thread)",
        #     "//",
        # ),
        # Bar(
        #     "dom_nl_24_20220220.csv",
        #     "gtc:gt:cpu_kfirst",
        #     "cyan",
        #     "GT4Py CPU (24 threads)",
        # ),
        # all gt4py backends
        Bar(
            "dom_nl_24_20220220.csv",
            "gtc:gt:cpu_kfirst",
            "cyan",
            "gtc:gt:cpu_kfirst",
        ),
        Bar(
            "dom_nl_24_20220220.csv",
            "gtc:gt:cpu_ifirst",
            "orange",
            "gtc:gt:cpu_ifirst",
        ),
        Bar(
            "dom_nl_24_20220220.csv",
            "gtc:dace",
            "red",
            "gtc:dace",
        ),
        Bar("dom_nl_24_20220220.csv", "gtc:gt:gpu", "green", "gtc:gt:gpu"),
        Bar("dom_nl_24_20220220.csv", "gtc:cuda", "purple", "gtc:cuda"),
    ],
)
pool2 = Pool(
    "Tangent-linear",
    [
        # fortran vs gt4py
        # Bar("dom_tl_24_20220222.csv", "fortran", "blue", "Fortran"),
        # Bar(
        #     "dom_tl_24_20220222.csv",
        #     "gtc:gt:cpu_kfirst",
        #     "cyan",
        #     "GT4Py CPU",
        # ),
        # Bar("dom_tl_24_20220222.csv", "gtc:gt:gpu", "green", "GT4Py GPU"),
        # all gt4py backends
        Bar(
            "dom_tl_24_20220222.csv",
            "gtc:gt:cpu_kfirst",
            "cyan",
            "gtc:gt:cpu_kfirst",
        ),
        Bar(
            "dom_tl_24_20220222.csv",
            "gtc:gt:cpu_ifirst",
            "orange",
            "gtc:gt:cpu_ifirst",
        ),
        Bar(
            "dom_tl_24_20220222.csv",
            "gtc:dace",
            "red",
            "gtc:dace",
        ),
        Bar("dom_tl_24_20220222.csv", "gtc:gt:gpu", "green", "gtc:gt:gpu"),
        Bar("dom_tl_24_20220222.csv", "gtc:cuda", "purple", "gtc:cuda"),
    ],
)
pool3 = Pool(
    "Adjoint",
    [
        # fortran vs gt4py
        # Bar("dom_ad_24_20220222.csv", "fortran", "blue", "Fortran"),
        # Bar(
        #     "dom_ad_24_20220222.csv",
        #     "gtc:gt:cpu_kfirst",
        #     "cyan",
        #     "GT4Py CPU",
        # ),
        # all gt4py backends
        Bar(
            "dom_ad_24_20220222.csv",
            "gtc:gt:cpu_kfirst",
            "cyan",
            "gtc:gt:cpu_kfirst",
        ),
        Bar(
            "dom_ad_24_20220222.csv",
            "gtc:gt:cpu_ifirst",
            "orange",
            "gtc:gt:cpu_ifirst",
        ),
        Bar(
            "dom_ad_24_20220222.csv",
            "gtc:dace",
            "red",
            "gtc:dace",
        ),
    ],
)
pools = [pool1, pool2, pool3]


figure_properties = {"figsize": [12, 6], "fontsize": 16, "tight_layout": True}
axes_properties = {
    "fontsize": 16,
    "x_lim": None,
    "x_ticks": None,
    "x_tick_length": 0,
    "y_label": "Runtime [ms]",
    "y_lim": [0, 1600],
    "y_ticks": None,
    "y_tick_format": "%4.d",
    "legend_on": True,
    "legend_fontsize": 14,
    "legend_framealpha": 1.0,
    "legend_loc": "upper left",
    "legend_ncol": 1,
}


def get_xs():
    axes_properties["x_ticks"] = []
    axes_properties["x_tick_labels"] = []
    xs = []
    start = 1
    for pool in pools:
        stop = start + len(pool.bars)
        xs.append(np.arange(start, stop))
        axes_properties["x_ticks"].append(0.5 * (start + stop - 1))
        axes_properties["x_tick_labels"].append(pool.name)
        start = stop + 0.75
    axes_properties["x_lim"] = [-0.25, start - 0.5]
    return xs


def plot_bars(ax, xs):
    for id_pool, pool in enumerate(pools):
        for bar, bar_center in zip(pool.bars, xs[id_pool]):
            df = pd.read_csv(os.path.join("timings", bar.csv_file))
            data = df[df[bar.column].notna()][bar.column]
            bar_height = sum(data) / len(data)
            ax.bar(
                bar_center,
                bar_height,
                width=1.0,
                color=bar.color,
                edgecolor="black",
                label=bar.label if id_pool == 0 else None,
                hatch=bar.hatch,
            )
            ax.annotate(
                f"{bar_height:.2f}",
                xy=(bar_center, bar_height),
                ha="center",
                va="bottom",
                fontsize=12,
            )


def main():
    fig, ax = plot_utils.get_figure_and_axes(**figure_properties)
    xs = get_xs()
    plot_bars(ax, xs)
    plot_utils.set_axes_properties(ax, **axes_properties)
    plot_utils.set_figure_properties(fig, **figure_properties)
    plt.show()


if __name__ == "__main__":
    main()
