#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MaxNLocator


# --------------------------------------------------
# User settings
# --------------------------------------------------
directory = "./"          # directory containing run_i folders
train_file = "train/learning-curve.out"

xlim = (0, 80)
energy_ylim = (0, 0.001)
force_ylim = (0, 0.04)

save_figures = True
show_figures = False
output_dir = "learning_curves"


# --------------------------------------------------
# Prepare folders
# --------------------------------------------------
folders = [
    f for f in os.listdir(directory)
    if os.path.isdir(os.path.join(directory, f))
]

folders = sorted(folders)

if save_figures:
    os.makedirs(output_dir, exist_ok=True)


# --------------------------------------------------
# Plot learning curves
# --------------------------------------------------
for folder in folders:
    filename = os.path.join(directory, folder, train_file)

    if not os.path.isfile(filename):
        print(f"Skipping {folder}: no {train_file}")
        continue

    print(f"Plotting {folder}")

    data = np.loadtxt(filename)

    epoch = data[:, 0]

    energy_train = data[:, 1]
    energy_test = data[:, 2]

    force_train = data[:, 9]
    force_test = data[:, 10]

    fig, ax_E = plt.subplots(figsize=(7, 4))

    ax_E.set_xlabel("Epoch")
    ax_E.set_ylabel("Energy RMSE")
    ax_E.set_xlim(*xlim)
    ax_E.set_ylim(*energy_ylim)
    ax_E.grid(True)

    ax_F = ax_E.twinx()
    ax_F.set_ylabel("Force RMSE")
    ax_F.set_ylim(*force_ylim)
    ax_F.xaxis.set_major_locator(MaxNLocator(integer=True))

    cmap = matplotlib.colormaps["RdYlGn"]

    lns1 = ax_E.plot(epoch, energy_train, "-", marker="o",
                     label="E Train", color=cmap(0.1))
    lns2 = ax_E.plot(epoch, energy_test, "-", marker="o",
                     label="E Test", color=cmap(0.4))

    lns3 = ax_F.plot(epoch, force_train, "--", marker="s",
                     label="F Train", color=cmap(0.7))
    lns4 = ax_F.plot(epoch, force_test, "--", marker="s",
                     label="F Test", color=cmap(0.9))

    lines = lns1 + lns2 + lns3 + lns4
    labels = [line.get_label() for line in lines]

    ax_E.legend(lines, labels)
    ax_E.set_title(folder)

    fig.tight_layout()

    if save_figures:
        outfile = os.path.join(output_dir, f"{folder}_learning_curve.png")
        fig.savefig(outfile, dpi=300)

    if show_figures:
        plt.show()

    plt.close(fig)
