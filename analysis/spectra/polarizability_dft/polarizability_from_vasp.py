#!/usr/bin/env python3
import numpy as np

# --- constants ---
EPS0 = 8.8541878128e-12        # vacuum permittivity (F/m)
ANGSTROM_TO_M = 1e-10
AU_POL = 1.64877727436e-41      # 1 atomic unit of polarizability in F·m²

# --- read dielectric tensor ---
eps_tensor = np.loadtxt("dial_tensor")
if eps_tensor.shape != (3, 3):
    raise ValueError("Dielectric tensor file must be 3x3 numbers.")

# --- read POSCAR lattice vectors ---
with open("POSCAR", "r") as f:
    lines = f.readlines()

scale = float(lines[1].strip())
a = np.array(list(map(float, lines[2].split())))
b = np.array(list(map(float, lines[3].split())))
c = np.array(list(map(float, lines[4].split())))

# cell volume in Å³ and m³
volume_A3 = abs(np.dot(a, np.cross(b, c))) * scale**3
volume_m3 = volume_A3 * (ANGSTROM_TO_M**3)

# --- compute polarizability tensor ---
alpha_tensor_SI = EPS0 * volume_m3 * (eps_tensor - np.eye(3))  # in F·m²
alpha_tensor_au = alpha_tensor_SI / AU_POL                     # in atomic units
alpha_iso = np.trace(alpha_tensor_au) / 3.0                    # isotropic average

# --- save to file ---
output_file = "polarizability_tensor.txt"
np.savetxt(
    output_file,
    alpha_tensor_au,
    fmt="%15.8f",
    header=f"Polarizability tensor (atomic units)\nCell volume: {volume_A3:.6f} Å³\nIsotropic mean: {alpha_iso:.6f} a.u."
)

print(f"Polarizability tensor written to '{output_file}'")
print("Isotropic mean α_iso = {:.6f} a.u.".format(alpha_iso))
