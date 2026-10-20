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

output_dir = "comparison_plots"
save_figures = True
show_figures = False


# =========================================================
# Helper functions
# =========================================================

def get_directories(path):
    folders = [
        entry.name for entry in os.scandir(path)
        if entry.is_dir()
    ]
    return sorted(folders)


def get_square_limits(x, y, margin=0.05):
    """Return symmetric square limits covering x and y."""
    vmin = min(np.min(x), np.min(y))
    vmax = max(np.max(x), np.max(y))

    delta = vmax - vmin
    if delta == 0:
        delta = abs(vmax) if vmax != 0 else 1.0

    pad = margin * delta

    return vmin - pad, vmax + pad


def parity_plot(x, y, title, xlabel, ylabel, outfile=None):
    """Reference-vs-prediction square parity plot."""

    lim_min, lim_max = get_square_limits(x, y)

    fig, ax = plt.subplots(figsize=(5, 5))

    ax.scatter(x, y, s=2, alpha=0.6)

    ax.plot([lim_min, lim_max], [lim_min, lim_max], "--", label="y = x")

    ax.set_xlim(lim_min, lim_max)
    ax.set_ylim(lim_min, lim_max)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True)
    ax.legend()

    fig.tight_layout()

    if outfile is not None:
        fig.savefig(outfile, dpi=300)

    if show_figures:
        plt.show()

    plt.close(fig)


def rmse(x, y):
    return np.sqrt(np.mean((x - y) ** 2))


# =========================================================
# Main analysis
# =========================================================

folders = get_directories(path)
print("Folders found:")
print(folders)

if save_figures:
    os.makedirs(output_dir, exist_ok=True)


for folder in folders:

    force_path = os.path.join(path, folder, forces_file)
    energy_path = os.path.join(path, folder, energy_file)

    if not os.path.isfile(force_path):
        print(f"Skipping {folder}: missing {forces_file}")
        continue

    if not os.path.isfile(energy_path):
        print(f"Skipping {folder}: missing {energy_file}")
        continue

    print(f"\nAnalyzing {folder}")

    forces = np.loadtxt(force_path).T
    energies = np.loadtxt(energy_path).T

    # -----------------------------------------------------
    # forces.comp columns
    # 0: structure index
    # 1: atom index
    # 2: reference force
    # 3: NN force
    # -----------------------------------------------------

    f_ref = forces[2]
    f_nn = forces[3]

    # -----------------------------------------------------
    # energy.comp columns
    # 0: structure index
    # 1: number of atoms
    # 2: reference total energy
    # 3: NN total energy
    #
    # Convert total energy to energy per atom
    # -----------------------------------------------------

    n_atoms = energies[1]
    e_ref = energies[2] / n_atoms
    e_nn = energies[3] / n_atoms

    e_rmse = rmse(e_ref, e_nn)
    f_rmse = rmse(f_ref, f_nn)

   

    # -----------------------------------------------------
# Combined E + F plot for this run
# -----------------------------------------------------

combined_outfile = None

if save_figures:
    combined_outfile = os.path.join(output_dir, f"{folder}_combined_parity.png")

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# Energy plot
lim_min, lim_max = get_square_limits(e_ref, e_nn)

axes[0].scatter(e_ref, e_nn, s=2, alpha=0.6)
axes[0].plot([lim_min, lim_max], [lim_min, lim_max], "--")
axes[0].set_xlim(lim_min, lim_max)
axes[0].set_ylim(lim_min, lim_max)
axes[0].set_aspect("equal", adjustable="box")
axes[0].set_xlabel("Reference energy / atom")
axes[0].set_ylabel("NN energy / atom")
axes[0].set_title("Energy")
axes[0].grid(True)

# Force plot
lim_min, lim_max = get_square_limits(f_ref, f_nn)

axes[1].scatter(f_ref, f_nn, s=2, alpha=0.6)
axes[1].plot([lim_min, lim_max], [lim_min, lim_max], "--")
axes[1].set_xlim(lim_min, lim_max)
axes[1].set_ylim(lim_min, lim_max)
axes[1].set_aspect("equal", adjustable="box")
axes[1].set_xlabel("Reference force")
axes[1].set_ylabel("NN force")
axes[1].set_title("Force")
axes[1].grid(True)

fig.suptitle(folder)
fig.tight_layout()

if save_figures:
    fig.savefig(combined_outfile, dpi=300)

if show_figures:
    plt.show()

plt.close(fig)

   print(f"Energy RMSE: {e_rmse:.8f} eV/atom")
    print(f"Force  RMSE: {f_rmse:.8f} eV/Angstrom")
    print(f"Energy RMSE: {e_rmse * 1000:.4f} meV/atom")
    print(f"Force  RMSE: {f_rmse * 1000:.4f} meV/Angstrom")

    # -----------------------------------------------------
    # Plot energies
    # -----------------------------------------------------

    energy_outfile = None
    force_outfile = None

    if save_figures:
        energy_outfile = os.path.join(output_dir, f"{folder}_energy_parity.png")
        force_outfile = os.path.join(output_dir, f"{folder}_force_parity.png")

    parity_plot(
        e_ref,
        e_nn,
        title=f"Energies {folder}",
        xlabel="Reference energy / atom",
        ylabel="NN energy / atom",
        outfile=energy_outfile,
    )

    parity_plot(
        f_ref,
        f_nn,
        title=f"Forces {folder}",
        xlabel="Reference force",
        ylabel="NN force",
        outfile=force_outfile,
    )
