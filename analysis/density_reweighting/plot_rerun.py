#!/usr/bin/env python3

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ==========================
# Settings
# ==========================

ROOT = Path(".")
NATOMS = None   # set e.g. 64, 128, 256 if you want energy per atom

FILES = {
    "E_mW": "PotEng_mW.data",
    "E_NNP": "PotEng_NNP.data",
    "P_mW": "Pvir_mW.data",
    "P_NNP": "Pvir_NNP.data",
}


# ==========================
# Helpers
# ==========================

def parse_density(folder_name):
    """
    Extract density from folder names like rho_0.050.
    """
    m = re.search(r"rho_([0-9.]+)", folder_name)
    if m is None:
        return None
    return float(m.group(1))


def load_array(path):
    arr = np.loadtxt(path)
    return np.atleast_1d(arr)


def metrics(true, pred):
    diff = pred - true
    return {
        "bias": np.mean(diff),
        "mae": np.mean(np.abs(diff)),
        "rmse": np.sqrt(np.mean(diff**2)),
        "std_diff": np.std(diff),
    }


# ==========================
# Find folders
# ==========================

rho_dirs = sorted(
    [p for p in ROOT.rglob("rho_*") if p.is_dir()],
    key=lambda p: parse_density(p.name) if parse_density(p.name) is not None else 999,
)

if not rho_dirs:
    raise RuntimeError("No rho_* folders found.")

rows = []
all_data = []

for d in rho_dirs:
    rho = parse_density(d.name)

    if rho is None:
        print(f"Skipping {d}: cannot parse density")
        continue

    missing = [fname for fname in FILES.values() if not (d / fname).exists()]
    if missing:
        print(f"Skipping {d}: missing {missing}")
        continue

    E_mW = load_array(d / FILES["E_mW"])
    E_NNP = load_array(d / FILES["E_NNP"])
    P_mW = load_array(d / FILES["P_mW"])
    P_NNP = load_array(d / FILES["P_NNP"])

    n = min(len(E_mW), len(E_NNP), len(P_mW), len(P_NNP))

    E_mW = E_mW[:n]
    E_NNP = E_NNP[:n]
    P_mW = P_mW[:n]
    P_NNP = P_NNP[:n]

    if NATOMS is not None:
        E_mW_plot = E_mW / NATOMS
        E_NNP_plot = E_NNP / NATOMS
    else:
        E_mW_plot = E_mW
        E_NNP_plot = E_NNP

    Pmet = metrics(P_mW, P_NNP)
    Emet = metrics(E_mW_plot, E_NNP_plot)

    rows.append({
        "folder": str(d),
        "rho": rho,
        "nframes": n,

        "P_mW_mean": np.mean(P_mW),
        "P_mW_std": np.std(P_mW),
        "P_NNP_mean": np.mean(P_NNP),
        "P_NNP_std": np.std(P_NNP),
        "P_bias": Pmet["bias"],
        "P_mae": Pmet["mae"],
        "P_rmse": Pmet["rmse"],
        "P_std_diff": Pmet["std_diff"],

        "E_mW_mean": np.mean(E_mW_plot),
        "E_mW_std": np.std(E_mW_plot),
        "E_NNP_mean": np.mean(E_NNP_plot),
        "E_NNP_std": np.std(E_NNP_plot),
        "E_bias": Emet["bias"],
        "E_mae": Emet["mae"],
        "E_rmse": Emet["rmse"],
    })

    for i in range(n):
        all_data.append({
            "rho": rho,
            "E_mW": E_mW_plot[i],
            "E_NNP": E_NNP_plot[i],
            "P_mW": P_mW[i],
            "P_NNP": P_NNP[i],
        })


summary = pd.DataFrame(rows).sort_values("rho")
data = pd.DataFrame(all_data)

summary.to_csv("density_scan_summary.csv", index=False)
data.to_csv("density_scan_all_frames.csv", index=False)

print(summary[[
    "rho", "nframes",
    "P_mW_mean", "P_NNP_mean", "P_bias", "P_rmse",
    "E_bias", "E_rmse"
]])


# ==========================
# Plot 1: EOS + bias summary
# ==========================

fig, axes = plt.subplots(2, 2, figsize=(11, 8))

# Pressure EOS
axes[0, 0].errorbar(
    summary["rho"], summary["P_mW_mean"],
    yerr=summary["P_mW_std"],
    marker="o", capsize=3, label="mW"
)
axes[0, 0].errorbar(
    summary["rho"], summary["P_NNP_mean"],
    yerr=summary["P_NNP_std"],
    marker="o", capsize=3, label="NNP"
)
axes[0, 0].axhline(0, linestyle="--", linewidth=1)
axes[0, 0].set_xlabel(r"Density $\rho$ (g/cm$^3$)")
axes[0, 0].set_ylabel("Virial pressure (bar)")
axes[0, 0].set_title("Pressure-density relation")
axes[0, 0].legend()
axes[0, 0].grid(alpha=0.3)

# Pressure bias
axes[0, 1].plot(summary["rho"], summary["P_bias"], marker="o")
axes[0, 1].axhline(0, linestyle="--", linewidth=1)
axes[0, 1].set_xlabel(r"Density $\rho$ (g/cm$^3$)")
axes[0, 1].set_ylabel(r"$P_{\rm NNP} - P_{\rm mW}$ (bar)")
axes[0, 1].set_title("Pressure bias")
axes[0, 1].grid(alpha=0.3)

# Pressure RMSE
axes[1, 0].plot(summary["rho"], summary["P_rmse"], marker="o")
axes[1, 0].set_xlabel(r"Density $\rho$ (g/cm$^3$)")
axes[1, 0].set_ylabel("Pressure RMSE (bar)")
axes[1, 0].set_title("Pressure error")
axes[1, 0].grid(alpha=0.3)

# Energy bias
axes[1, 1].plot(summary["rho"], summary["E_bias"], marker="o", label="bias")
axes[1, 1].plot(summary["rho"], summary["E_rmse"], marker="o", label="RMSE")
axes[1, 1].axhline(0, linestyle="--", linewidth=1)
axes[1, 1].set_xlabel(r"Density $\rho$ (g/cm$^3$)")
if NATOMS is None:
    axes[1, 1].set_ylabel("Energy error (eV)")
else:
    axes[1, 1].set_ylabel("Energy error per atom (eV)")
axes[1, 1].set_title("Energy error")
axes[1, 1].legend()
axes[1, 1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig("density_scan_summary.png", dpi=300)
plt.show()


# ==========================
# Plot 2: pressure parity by density
# ==========================

fig, ax = plt.subplots(figsize=(6, 6))

for rho, group in data.groupby("rho"):
    ax.scatter(
        group["P_mW"], group["P_NNP"],
        s=14, alpha=0.55, label=f"{rho:.3f}"
    )

lo = min(data["P_mW"].min(), data["P_NNP"].min())
hi = max(data["P_mW"].max(), data["P_NNP"].max())
pad = 0.05 * (hi - lo)

ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], linestyle="--", linewidth=1)
ax.set_xlim(lo - pad, hi + pad)
ax.set_ylim(lo - pad, hi + pad)

ax.set_xlabel("mW virial pressure (bar)")
ax.set_ylabel("NNP virial pressure (bar)")
ax.set_title("Pressure parity by density")
ax.legend(title=r"$\rho$", fontsize=8)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("pressure_parity_by_density.png", dpi=300)
plt.show()
