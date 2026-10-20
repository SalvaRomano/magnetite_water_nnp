#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt


# =========================================================
# User settings
# =========================================================

path = "./"

forces_file = "train/dataset/forces.comp"
energy_file = "train/dataset/energy.comp"

output_file = "all_energy_force_parity.png"

save_figure = True
show_figure = True


# =========================================================
# Helper functions
# =========================================================

def get_directories(path):
    return sorted([
        entry.name for entry in os.scandir(path)
        if entry.is_dir()
    ])


def get_square_limits(x, y, margin=0.05):
    vmin = min(np.min(x), np.min(y))
    vmax = max(np.max(x), np.max(y))

    delta = vmax - vmin
    if delta == 0:
        delta = abs(vmax) if vmax != 0 else 1.0

    pad = margin * delta
    return vmin - pad, vmax + pad


def rmse(x, y):
    return np.sqrt(np.mean((x - y) ** 2))


def parity_axis(ax, x, y, title, xlabel, ylabel, size=2):
    lim_min, lim_max = get_square_limits(x, y)

    ax.scatter(x, y, s=size, alpha=0.5)
    ax.plot([lim_min, lim_max], [lim_min, lim_max], "--")

    ax.set_xlim(lim_min, lim_max)
    ax.set_ylim(lim_min, lim_max)
    ax.set_aspect("equal", adjustable="box")

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True)


# =========================================================
# Load data
# =========================================================

folders = get_directories(path)

runs = []

for folder in folders:

    force_path = os.path.join(path, folder, forces_file)
    energy_path = os.path.join(path, folder, energy_file)

    if not os.path.isfile(force_path) or not os.path.isfile(energy_path):
        print(f"Skipping {folder}: missing force or energy file")
        continue

    forces = np.loadtxt(force_path).T
    energies = np.loadtxt(energy_path).T

    f_ref = forces[2]
    f_nn = forces[3]

    n_atoms = energies[1]
    e_ref = energies[2] / n_atoms
    e_nn = energies[3] / n_atoms

    runs.append({
        "folder": folder,
        "e_ref": e_ref,
        "e_nn": e_nn,
        "f_ref": f_ref,
        "f_nn": f_nn,
        "e_rmse": rmse(e_ref, e_nn),
        "f_rmse": rmse(f_ref, f_nn),
    })


# =========================================================
# Plot all runs in one large figure
# =========================================================

n_runs = len(runs)

if n_runs == 0:
    raise RuntimeError("No valid runs found.")

fig, axes = plt.subplots(
    n_runs,
    2,
    figsize=(10, 4 * n_runs),
    squeeze=False
)

for i, run in enumerate(runs):

    folder = run["folder"]

    e_title = (
        f"{folder} | Energy "
        f"RMSE = {run['e_rmse'] * 1000:.3f} meV/atom"
    )

    f_title = (
        f"{folder} | Force "
        f"RMSE = {run['f_rmse'] * 1000:.3f} meV/Å"
    )

    parity_axis(
        axes[i, 0],
        run["e_ref"],
        run["e_nn"],
        title=e_title,
        xlabel="Reference energy / atom",
        ylabel="NN energy / atom",
        size=3,
    )

    parity_axis(
        axes[i, 1],
        run["f_ref"],
        run["f_nn"],
        title=f_title,
        xlabel="Reference force",
        ylabel="NN force",
        size=1,
    )

fig.tight_layout()

if save_figure:
    fig.savefig(output_file, dpi=300)

if show_figure:
    plt.show()
else:
    plt.close(fig)
