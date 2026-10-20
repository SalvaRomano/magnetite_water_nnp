#!/usr/bin/env python3

def parse_lammpstrj(traj_file):
    """
    Parse a LAMMPS trajectory and yield one frame at a time.
    Assumes a block structure like:

    ITEM: TIMESTEP
    ...
    ITEM: NUMBER OF ATOMS
    ...
    ITEM: BOX BOUNDS pp pp pp
    xlo xhi
    ylo yhi
    zlo zhi
    ITEM: ATOMS id type x y z fx fy fz
    ...
    """
    with open(traj_file, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break

            line = line.strip()
            if line != "ITEM: TIMESTEP":
                continue

            timestep = int(f.readline().strip())

            line = f.readline().strip()
            if line != "ITEM: NUMBER OF ATOMS":
                raise ValueError(f"Expected 'ITEM: NUMBER OF ATOMS', got: {line}")
            natoms = int(f.readline().strip())

            line = f.readline().strip()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"Expected 'ITEM: BOX BOUNDS ...', got: {line}")

            xlo, xhi = map(float, f.readline().split())
            ylo, yhi = map(float, f.readline().split())
            zlo, zhi = map(float, f.readline().split())

            lx = xhi - xlo
            ly = yhi - ylo
            lz = zhi - zlo

            line = f.readline().strip()
            if not line.startswith("ITEM: ATOMS"):
                raise ValueError(f"Expected 'ITEM: ATOMS ...', got: {line}")

            columns = line.split()[2:]
            col_index = {name: i for i, name in enumerate(columns)}

            required = ["x", "y", "z", "fx", "fy", "fz"]
            missing = [c for c in required if c not in col_index]
            if missing:
                raise ValueError(f"Missing required columns in trajectory: {missing}")

            atoms = []
            for _ in range(natoms):
                parts = f.readline().split()
                x = float(parts[col_index["x"]])
                y = float(parts[col_index["y"]])
                z = float(parts[col_index["z"]])
                fx = float(parts[col_index["fx"]])
                fy = float(parts[col_index["fy"]])
                fz = float(parts[col_index["fz"]])
                atoms.append((x, y, z, fx, fy, fz))

            yield {
                "timestep": timestep,
                "natoms": natoms,
                "xlo": xlo,
                "ylo": ylo,
                "zlo": zlo,
                "lx": lx,
                "ly": ly,
                "lz": lz,
                "atoms": atoms,
            }


def read_energies(energy_file):
    """Read one energy per non-empty line."""
    energies = []
    with open(energy_file, "r") as f:
        for line in f:
            s = line.strip()
            if s:
                energies.append(float(s))
    return energies


def write_input_data(frames, energies, output_file, atom_label="O"):
    """
    Write one block per frame in the format:

    begin
    comment
    lattice ...
    lattice ...
    lattice ...
    atom x y z O <c1> <n1> fx fy fz
    ...
    energy ...
    charge <unused>
    end
    """
    frames = list(frames)

    if len(frames) != len(energies):
        raise ValueError(
            f"Number of frames ({len(frames)}) != number of energies ({len(energies)})"
        )

    with open(output_file, "w") as out:
        for frame, energy in zip(frames, energies):
            out.write("begin\n")
            out.write("comment\n")

            out.write(f"lattice {frame['lx']:.14f} 0.000000000000000E+000 0.000000000000000E+000\n")
            out.write(f"lattice 0.000000000000000E+000 {frame['ly']:.14f} 0.000000000000000E+000\n")
            out.write(f"lattice 0.000000000000000E+000 0.000000000000000E+000 {frame['lz']:.14f}\n")

            xlo = frame["xlo"]
            ylo = frame["ylo"]
            zlo = frame["zlo"]

            for x, y, z, fx, fy, fz in frame["atoms"]:
                x_shift = x - xlo
                y_shift = y - ylo
                z_shift = z - zlo

                out.write(
                    f"atom {x_shift:.14f} {y_shift:.14f} {z_shift:.14f} "
                    f"{atom_label} <c1> <n1> "
                    f"{fx:.15f} {fy:.15f} {fz:.15f}\n"
                )

            out.write(f"energy {energy:.15f}\n")
            out.write("charge <unused>\n")
            out.write("end\n")


def main():
    traj_file = "water.lammpstrj"
    energy_file = "E.data"
    output_file = "input.data"

    frames = parse_lammpstrj(traj_file)
    energies = read_energies(energy_file)
    write_input_data(frames, energies, output_file)

    print(f"Wrote {output_file}")


if __name__ == "__main__":
    main()
