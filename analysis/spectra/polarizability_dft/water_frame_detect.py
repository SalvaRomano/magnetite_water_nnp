#!/usr/bin/env python3
import numpy as np
import sys

def normalize(v):
    v = np.array(v, dtype=float)
    n = np.linalg.norm(v)
    return v / n if n != 0 else v

def read_poscar(filename):
    with open(filename) as f:
        lines = [l.strip() for l in f if l.strip()]
    scale = float(lines[1])
    lattice = np.array([list(map(float, lines[i].split())) for i in range(2,5)]) * scale
    species = lines[5].split()
    counts = list(map(int, lines[6].split()))
    n_atoms = sum(counts)
    coord_type = lines[7].lower()
    coords = np.array([list(map(float, lines[8+i].split()[:3])) for i in range(n_atoms)])
    if coord_type.startswith('d'):  # Direct
        coords = coords @ lattice
    elif coord_type.startswith('c'):
        coords *= scale
    else:
        raise ValueError("Unknown coordinate type (expect Direct or Cartesian)")
    return species, counts, coords, lattice

def identify_water(coords, species, counts, H_index, cutoff=1.2):
    atom_types = []
    for sp, n in zip(species, counts):
        atom_types += [sp]*n
    if atom_types[H_index].upper() != 'H':
        raise ValueError("Selected atom is not hydrogen!")
    H_pos = coords[H_index]
    oxy_indices = [i for i,t in enumerate(atom_types) if t.upper() == 'O']
    hyd_indices = [i for i,t in enumerate(atom_types) if t.upper() == 'H']
    dists = np.linalg.norm(coords[oxy_indices] - H_pos, axis=1)
    O_index = oxy_indices[np.argmin(dists)]
    O_pos = coords[O_index]
    if dists.min() > cutoff:
        raise ValueError("No oxygen found within reasonable O–H distance.")
    dH = np.linalg.norm(coords[hyd_indices] - O_pos, axis=1)
    near_H = [hyd_indices[i] for i,d in enumerate(dH) if d < cutoff]
    if len(near_H) != 2:
        raise ValueError("Expected exactly two hydrogens near this O.")
    H_other_index = near_H[0] if near_H[1]==H_index else near_H[1]
    H_other = coords[H_other_index]
    return O_index, H_other_index, O_pos, H_pos, H_other

def rotation_matrix(O, H_sel, H_other):
    z_axis = normalize(H_sel - O)
    bis = normalize((H_sel - O) + (H_other - O))
    x_axis = normalize(bis)
    y_axis = normalize(np.cross(z_axis, x_axis))
    x_axis = normalize(np.cross(y_axis, z_axis))
    R = np.vstack([x_axis, y_axis, z_axis])
    return R

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 water_frame_detect.py POSCAR_ref POSCAR_disp")
        sys.exit(1)

    poscar_ref, poscar_disp = sys.argv[1], sys.argv[2]

    # --- read both POSCAR files ---
    species, counts, coords_ref, _ = read_poscar(poscar_ref)
    _, _, coords_disp, _ = read_poscar(poscar_disp)

    atom_types = []
    for sp, n in zip(species, counts):
        atom_types += [sp]*n

    # --- identify displaced hydrogen ---
    disp = np.linalg.norm(coords_disp - coords_ref, axis=1)
    H_indices = [i for i,t in enumerate(atom_types) if t.upper() == 'H']
    disp_H = disp[H_indices]
    displaced_H_index = H_indices[np.argmax(disp_H)]
    print(f"Detected displaced hydrogen index: {displaced_H_index} (0-based)")

    # --- build molecular frame based on that H ---
    O_idx, H_other_idx, O_pos, H_sel, H_other = identify_water(coords_ref, species, counts, displaced_H_index)
    R = rotation_matrix(O_pos, H_sel, H_other)

    # --- save rotation matrix ---
    header = (
        f"Rotation matrix lab→molecular frame (rows = x,y,z)\n"
        f"Displaced H index: {displaced_H_index}\n"
        f"O index: {O_idx}, other H index: {H_other_idx}"
    )
    np.savetxt("rotation_matrix.txt", R, fmt="%12.8f", header=header)
    print("Rotation matrix written to 'rotation_matrix.txt'.")

if __name__ == "__main__":
    main()
