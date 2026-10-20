#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

AVOGADRO = 6.02214076e23  # mol^-1
ANG3_TO_CM3 = 1e-24


def read_lammps_data(filename):
    filename = Path(filename)
    lines = filename.read_text().splitlines()

    # ---- number of atoms ----
    natoms = None
    for line in lines:
        m = re.match(r"\s*(\d+)\s+atoms", line)
        if m:
            natoms = int(m.group(1))
            break

    if natoms is None:
        raise ValueError(f"Could not find number of atoms in {filename}")

    # ---- box bounds ----
    bounds = {}
    for i, line in enumerate(lines):
        if "xlo xhi" in line:
            xlo, xhi = map(float, line.split()[:2])
            bounds["x"] = [xlo, xhi, i]
        elif "ylo yhi" in line:
            ylo, yhi = map(float, line.split()[:2])
            bounds["y"] = [ylo, yhi, i]
        elif "zlo zhi" in line:
            zlo, zhi = map(float, line.split()[:2])
            bounds["z"] = [zlo, zhi, i]

    if len(bounds) != 3:
        raise ValueError(f"Could not find full box bounds in {filename}")

    # ---- masses ----
    masses = {}
    in_masses = False

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("Masses"):
            in_masses = True
            continue

        if in_masses:
            if stripped == "":
                continue

            parts = stripped.split()

            # stop if next section starts
            if not parts[0].isdigit():
                break

            atom_type = int(parts[0])
            mass = float(parts[1])
            masses[atom_type] = mass

    if not masses:
        raise ValueError(f"Could not read Masses section in {filename}")

    # ---- atoms ----
    atoms_start = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Atoms"):
            atoms_start = i
            break

    if atoms_start is None:
        raise ValueError(f"Could not find Atoms section in {filename}")

    # skip empty lines after "Atoms"
    first_atom_line = atoms_start + 1
    while lines[first_atom_line].strip() == "":
        first_atom_line += 1

    atom_lines = lines[first_atom_line:first_atom_line + natoms]

    atoms = []
    for line in atom_lines:
        parts = line.split()
        atom_id = int(parts[0])
        atom_type = int(parts[1])
        x, y, z = map(float, parts[2:5])

        atoms.append(
            {
                "id": atom_id,
                "type": atom_type,
                "x": x,
                "y": y,
                "z": z,
                "parts": parts,
            }
        )

    return lines, bounds, masses, atoms, first_atom_line, natoms


def total_mass_grams(atoms, masses):
    total_mass_amu = sum(masses[a["type"]] for a in atoms)
    return total_mass_amu / AVOGADRO


def box_volume_ang3(bounds):
    lx = bounds["x"][1] - bounds["x"][0]
    ly = bounds["y"][1] - bounds["y"][0]
    lz = bounds["z"][1] - bounds["z"][0]
    return lx * ly * lz


def density_g_cm3(atoms, masses, bounds):
    mass_g = total_mass_grams(atoms, masses)
    volume_cm3 = box_volume_ang3(bounds) * ANG3_TO_CM3
    return mass_g / volume_cm3


def write_scaled_data(input_file, output_file, target_rho):
    lines, bounds, masses, atoms, first_atom_line, natoms = read_lammps_data(input_file)

    rho0 = density_g_cm3(atoms, masses, bounds)

    scale = (rho0 / target_rho) ** (1.0 / 3.0)

    old_lx = bounds["x"][1] - bounds["x"][0]
    old_ly = bounds["y"][1] - bounds["y"][0]
    old_lz = bounds["z"][1] - bounds["z"][0]

    new_lx = old_lx * scale
    new_ly = old_ly * scale
    new_lz = old_lz * scale

    xlo, ylo, zlo = bounds["x"][0], bounds["y"][0], bounds["z"][0]

    new_bounds = {
        "x": [xlo, xlo + new_lx],
        "y": [ylo, ylo + new_ly],
        "z": [zlo, zlo + new_lz],
    }

    new_lines = lines.copy()

    # update box
    ix = bounds["x"][2]
    iy = bounds["y"][2]
    iz = bounds["z"][2]

    new_lines[ix] = f"{new_bounds['x'][0]:20.12f} {new_bounds['x'][1]:20.12f}  xlo xhi"
    new_lines[iy] = f"{new_bounds['y'][0]:20.12f} {new_bounds['y'][1]:20.12f}  ylo yhi"
    new_lines[iz] = f"{new_bounds['z'][0]:20.12f} {new_bounds['z'][1]:20.12f}  zlo zhi"

    # update coordinates
    for i, atom in enumerate(atoms):
        parts = atom["parts"]

        x_scaled = xlo + (atom["x"] - xlo) * scale
        y_scaled = ylo + (atom["y"] - ylo) * scale
        z_scaled = zlo + (atom["z"] - zlo) * scale

        parts[2] = f"{x_scaled:.12f}"
        parts[3] = f"{y_scaled:.12f}"
        parts[4] = f"{z_scaled:.12f}"

        new_lines[first_atom_line + i] = (
            f"{int(parts[0]):10d} {int(parts[1]):4d} "
            f"{float(parts[2]):20.12f} "
            f"{float(parts[3]):20.12f} "
            f"{float(parts[4]):20.12f}"
        )

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text("\n".join(new_lines) + "\n")

    return rho0, scale, new_lx, new_ly, new_lz


def main():
    parser = argparse.ArgumentParser(
        description="Create stretched LAMMPS data files at lower target densities."
    )

    parser.add_argument(
        "files",
        nargs="+",
        help="Input LAMMPS data files, e.g. N_64/initial.data N_128/initial.data",
    )

    parser.add_argument(
        "--rho",
        nargs="+",
        type=float,
        default=[0.8, 0.6, 0.4, 0.2, 0.1, 0.05],
        help="Target densities in g/cm^3.",
    )

    parser.add_argument(
        "--out",
        default="stretched_data",
        help="Output directory.",
    )

    args = parser.parse_args()

    for f in args.files:
        f = Path(f)
        system_name = f.parent.name if f.parent.name != "." else f.stem

        for rho in args.rho:
            out_file = Path(args.out) / system_name / f"rho_{rho:.3f}" / "initial.data"

            rho0, scale, lx, ly, lz = write_scaled_data(
                input_file=f,
                output_file=out_file,
                target_rho=rho,
            )

            print(
                f"{f} -> {out_file} | "
                f"rho0={rho0:.4f} g/cm^3, "
                f"rho={rho:.4f} g/cm^3, "
                f"scale={scale:.4f}, "
                f"L=({lx:.3f}, {ly:.3f}, {lz:.3f}) Å"
            )


if __name__ == "__main__":
    main()
