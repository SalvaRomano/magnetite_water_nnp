"""Create the final SFG spectra plot with legends and surface/bulk marker."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
import numpy as np


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "spectra_with_legend.pdf"

FILES = [
    "FT_2.dat",
    "FT_5.dat",
    "FT_6.dat",
    "FT_4.dat",
    "FT_3.dat",
    "FT_1_third.dat",
    "FT_1_second.dat",
    "FT_1.dat",
]
LABELS = [
    "oct-H2O",
    "tet1-H2O",
    "tet2-H2O",
    r"OH$^{-}$",
    "oct-OH",
    "third",
    "second",
    "Total",
]
LINESTYLES = ["-", "-", "-", "-", "-", "-.", "--", "-"]
COLORS = [
    "tab:green",
    "deepskyblue",
    "mediumorchid",
    "tab:blue",
    "tab:orange",
    "black",
    "black",
    "black",
]

X_MIN, X_MAX = 2000, 4000
Y_MIN, Y_MAX = -1.4, 1.0
Y_SCALE = 1e-5


def load_spectra():
    """Load spectra and apply the normalization used in the original script."""
    spectra = {}
    for filename, label in zip(FILES, LABELS):
        data = np.loadtxt(HERE / filename)
        data[:, 1] = -data[:, 1] / 150

        if filename in {"FT_4.dat", "FT_1.dat", "FT_1_third.dat", "FT_1_second.dat"}:
            data[:, 1] /= 2
        if filename == "FT_3.dat":
            data[:, 1] *= 2

        spectra[label] = data
    return spectra


def make_figure():
    plt.rcParams.update(plt.rcParamsDefault)
    cm_to_inch = 0.393701
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.linewidth": 0.7,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "lines.linewidth": 0.9,
            "mathtext.default": "regular",
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    spectra = load_spectra()
    # Wide source aspect; the manuscript scales this to the column width.
    fig, ax = plt.subplots(figsize=(9 * cm_to_inch, 5.2 * cm_to_inch), dpi=300)

    for label, linestyle, color in zip(LABELS, LINESTYLES, COLORS):
        data = spectra[label]
        ax.plot(data[:, 0], data[:, 1] * Y_SCALE, color=color, linestyle=linestyle)

    ax.set(
        xlim=(X_MIN, X_MAX),
        ylim=(Y_MIN, Y_MAX),
        xlabel=r"Wavenumber [cm$^{-1}$]",
        ylabel=r"Im$\chi^{(2)}$ [a.u.]",
    )
    ax.set_xticks([2000, 2500, 3000, 3500, 4000])
    ax.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    # One minor tick between adjacent major ticks, as in Figure 6.
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(axis="x", which="minor", length=2.0, width=0.5)
    ax.grid(True, color="0.68", linewidth=0.8)

    # Black curves: distance from the surface encoded by line style.
    distance_handles = [
        Line2D([], [], color="black", lw=1.1, linestyle="-", label=r"first 3.7 $\AA$"),
        Line2D([], [], color="black", lw=1.1, linestyle="--", label=r"second 7.2 $\AA$"),
        Line2D([], [], color="black", lw=1.1, linestyle="-.", label=r"third 10.6 $\AA$"),
    ]
    distance_legend = ax.legend(
        handles=distance_handles,
        loc="upper left",
        frameon=False,
        fontsize=7.5,
        handlelength=1.25,
        handletextpad=0.55,
        borderpad=0.15,
        labelspacing=0.25,
    )
    ax.add_artist(distance_legend)

    # Colored curves: chemical contributions. tet1 and tet2 share one entry.
    species_handles = [
        Line2D([], [], color="mediumorchid", lw=1.1, label=r"H$_2$O tet"),
        Line2D([], [], color="tab:green", lw=1.1, label=r"H$_2$O oct"),
        Line2D([], [], color="tab:orange", lw=1.1, label="remaining OH"),
        Line2D([], [], color="tab:blue", lw=1.1, label="dissociated OH"),
        Line2D([], [], color="deepskyblue", lw=1.1, label=r"H$_2$O on dissociation"),
    ]
    ax.legend(
        handles=species_handles,
        loc="lower left",
        frameon=False,
        fontsize=7.5,
        handlelength=1.25,
        handletextpad=0.55,
        borderpad=0.15,
        labelspacing=0.12,
    )

    # Direction guide reproduced from the final Inkscape figure.
    arrow = dict(arrowstyle="->", color="black", lw=0.8, shrinkA=0, shrinkB=0)
    ax.annotate("", xy=(1.035, 0.91), xytext=(1.035, 0.61),
                xycoords="axes fraction", arrowprops=arrow, annotation_clip=False)
    ax.annotate("", xy=(1.035, 0.20), xytext=(1.035, 0.54),
                xycoords="axes fraction", arrowprops=arrow, annotation_clip=False)
    ax.text(1.075, 0.76, "surface", rotation=90, va="center", ha="center",
            transform=ax.transAxes, fontsize=8)
    ax.text(1.075, 0.37, "bulk", rotation=90, va="center", ha="center",
            transform=ax.transAxes, fontsize=8)

    # Tight export removes the unused white frame while retaining annotations.
    fig.subplots_adjust(left=0.17, right=0.82, bottom=0.20, top=0.97)
    save_options = {"bbox_inches": "tight", "pad_inches": 0.01}
    fig.savefig(OUTPUT, **save_options)
    fig.savefig(OUTPUT.with_suffix(".svg"), **save_options)
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=600, **save_options)
    return fig


if __name__ == "__main__":
    make_figure()
    plt.show()
    print(OUTPUT)
