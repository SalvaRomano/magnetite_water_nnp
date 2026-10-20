#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

dt_ps = 0.001
out_freq = 5000   # change if needed

E = np.loadtxt("E.data")
L = np.loadtxt("L.data")
P = np.loadtxt("P.data")
rho = np.loadtxt("rho.data")
T = np.loadtxt("T.data")

n = len(E)
for name, arr in {"L": L, "P": P, "rho": rho, "T": T}.items():
    if len(arr) != n:
        raise ValueError(f"{name}.data has length {len(arr)}, expected {n}")

time_ns = np.arange(n) * out_freq * dt_ps / 1000.0

fig, axes = plt.subplots(2, 3, figsize=(10, 8), sharex=True)
axes = axes.ravel()

axes[0].plot(time_ns, E)
axes[0].set_ylabel("Epot (eV)")
axes[0].set_title("Potential energy")

axes[1].plot(time_ns, T)
axes[1].set_ylabel("T (K)")
axes[1].set_title("Temperature")

axes[2].plot(time_ns, P)
axes[2].set_ylabel("P (bar)")
axes[2].set_title("Pressure")

axes[3].plot(time_ns, rho)
axes[3].set_ylabel(r"$\rho$ (g/cm$^3$)")
axes[3].set_title("Density")

axes[4].plot(time_ns, L)
axes[4].set_ylabel(r"L ($\AA$)")
axes[4].set_title("Box length")

# leave last panel empty, or use it for text summary
axes[5].axis("off")
summary = (
    f"frames: {n}\n"
    f"mean T: {np.mean(T[:]):.1f} K\n"
    f"mean P: {np.mean(P[:]):.1f} bar\n"
    f"mean rho: {np.mean(rho[:]):.4f} g/cm$^3$\n"
    f"mean L: {np.mean(L[:]):.3f} Å"
)
axes[5].text(0.05, 0.7, summary, fontsize=11, va="top")

for ax in axes[:5]:
    ax.set_xlabel("Time (ns)")
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("npt_diagnostics_grid.png", dpi=300)
plt.show()
