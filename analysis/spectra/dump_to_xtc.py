#!/usr/bin/env python3

import argparse
import sys
import MDAnalysis as mda


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert LAMMPS dump (metal units) directly to "
            "traj.xtc (positions) and vel.xtc (velocities), "
            "auto-generating a minimal topology."
        )
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="LAMMPS dump file (must contain x y z vx vy vz)"
    )

    parser.add_argument(
        "--traj-out",
        default="traj.xtc",
        help="Output XTC for positions (default: traj.xtc)"
    )

    parser.add_argument(
        "--vel-out",
        default="vel.xtc",
        help="Output XTC for velocities (default: vel.xtc)"
    )

    parser.add_argument(
        "--top-out",
        default=None,
        help="Optional: write guessed topology to this PDB file"
    )

    args = parser.parse_args()

    # --- Load dump ---
    try:
        u = mda.Universe(args.input, format="LAMMPSDUMP")
    except Exception as e:
        sys.exit(f"ERROR loading dump file: {e}")

    n_atoms = u.atoms.n_atoms

    # --- Guess minimal topology ---
    # Required for writers and many tools
    u.add_TopologyAttr("name", ["X"] * n_atoms)
    u.add_TopologyAttr("type", ["X"] * n_atoms)
    u.add_TopologyAttr("resname", ["SYS"] * n_atoms)
    u.add_TopologyAttr("resid", [1] * n_atoms)

    # --- Optional: write topology for inspection ---
    if args.top_out is not None:
        u.trajectory[0]
        u.atoms.positions *= 0.1  # Å -> nm
        with mda.Writer(args.top_out) as W:
            W.write(u)
        print(f"Minimal topology written to: {args.top_out}")

    # --- Writers ---
    Wpos = mda.Writer(args.traj_out, n_atoms)
    Wvel = mda.Writer(args.vel_out, n_atoms, velocities=True)

    # --- Stream frames ---
    for ts in u.trajectory:
        # Convert units (metal → GROMACS)
        u.atoms.positions *= 0.1       # Å -> nm

        if ts.has_velocities:
            u.atoms.velocities *= 0.1  # Å/ps -> nm/ps
        else:
            sys.exit("ERROR: velocities not found in dump.")

        Wpos.write(u)
        Wvel.write(u)

    Wpos.close()
    Wvel.close()

    print("Conversion completed successfully:")
    print(f"  Positions  -> {args.traj_out}")
    print(f"  Velocities -> {args.vel_out}")
    print(f"  Atoms      -> {n_atoms}")


if __name__ == "__main__":
    main()

