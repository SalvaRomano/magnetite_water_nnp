#!/usr/bin/env python3

import os
import numpy as np


# =========================================================
# User settings
# =========================================================

directory = "./"
learning_file = "train/learning-curve.out"

# columns in learning-curve.out
epoch_col = 0

# Training errors
energy_col = 1
force_col = 9

# If you want test errors instead, use:
# energy_col = 2
# force_col = 10

w_energy = 0.3
w_force = 0.7

# Ignore very early epochs if needed
min_epoch = 0

# Optional: only consider epochs after this many epochs
max_epoch = None   # e.g. 80, or None

print("Best weights: to put in dataset.sh")


# =========================================================
# Collect folders
# =========================================================

folders = [
    f for f in os.listdir(directory)
    if os.path.isdir(os.path.join(directory, f))
]

folders = sorted(folders)


# =========================================================
# Select best epoch per run
# =========================================================

best_epochs = []
err_energy = []
err_force = []

for folder in folders:

    filename = os.path.join(directory, folder, learning_file)

    if not os.path.isfile(filename):
        print(f"\nSkipping {folder}: file not found")
        continue

    data = np.loadtxt(filename)

    epoch = data[:, epoch_col]
    energy = data[:, energy_col]
    force = data[:, force_col]

    mask = epoch >= min_epoch

    if max_epoch is not None:
        mask &= epoch <= max_epoch

    epoch_sel = epoch[mask]
    energy_sel = energy[mask]
    force_sel = force[mask]

    # Normalize by the best value of each error
    energy_norm = energy_sel / np.min(energy_sel)
    force_norm = force_sel / np.min(force_sel)

    score = w_energy * energy_norm + w_force * force_norm

    best_local = np.argmin(score)

    best_epoch = int(epoch_sel[best_local])
    best_E = energy_sel[best_local]
    best_F = force_sel[best_local]

    best_epochs.append(best_epoch)
    err_energy.append(best_E)
    err_force.append(best_F)

    print(best_epoch, end=" ")


# =========================================================
# Print summary
# =========================================================

print("\n")
print("Selected epochs:")
print(best_epochs)

print("\nE errors:")
print(np.round(np.array(err_energy), 8))

print("\nF errors:")
print(np.round(np.array(err_force), 8))

print("\nE errors (meV/atom):")
print(np.round(np.array(err_energy) * 1000, 5))

print("\nF errors (meV/Angstrom):")
print(np.round(np.array(err_force) * 1000, 4))
