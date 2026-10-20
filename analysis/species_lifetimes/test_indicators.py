#!/usr/bin/env python3
"""
test_indicators.py
-------------------
Read a (small) oxy_table.csv, extract species of interest,
build indicator functions, and print a preview.

Run:
    python test_indicators.py oxy_table.csv
"""

import sys
import pandas as pd
import numpy as np
from collections import defaultdict

# species we track for correlation
TARGET_SPECIES = ["dissociated_OH", "remaining_OH", "tet_closest"]


# =====================================================================
# Stage 1: Read CSV and build frame → set(atom_index) maps
# =====================================================================
def build_species_maps(csv_path: str):
    """
    Return:
        maps: {species: {frame: set(atom_indices)}}
        n_frames: total number of frames (max_frame + 1)
    """
    maps = {s: defaultdict(set) for s in TARGET_SPECIES}
    max_frame = -1

    usecols = ["frame", "atom_index", "species"]

    for chunk in pd.read_csv(csv_path, chunksize=200000, usecols=usecols):
        for _, row in chunk.iterrows():
            sp = row["species"]
            if sp in maps:
                f = int(row["frame"])
                a = int(row["atom_index"])
                maps[sp][f].add(a)
                if f > max_frame:
                    max_frame = f

    return maps, max_frame + 1


# =====================================================================
# Stage 2: Build indicator vectors for each oxygen that ever appears
# =====================================================================
def build_indicator_vectors(species_maps, n_frames):
    """
    indicators[species][atom_index] = boolean vector of length n_frames.
    """
    indicators = {s: {} for s in species_maps}

    for sp, frame_dict in species_maps.items():

        # collect all atoms ever observed in this species
        if frame_dict:
            all_atoms = set().union(*frame_dict.values())
        else:
            all_atoms = set()

        for atom in sorted(all_atoms):
            v = np.zeros(n_frames, dtype=np.uint8)
            for f, atoms in frame_dict.items():
                if atom in atoms:
                    v[f] = 1
            indicators[sp][atom] = v

    return indicators


# =====================================================================
# Pretty-Print Preview
# =====================================================================
def preview_indicators(indicators, max_atoms=15, max_frames=40):
    """
    Print a small readable preview of the indicator vectors.
    """
    print("\n====== INDICATOR FUNCTION PREVIEW ======\n")

    for sp, atom_map in indicators.items():
        print(f"Species: {sp}  (total atoms: {len(atom_map)})")

        atoms = list(atom_map.keys())
        if len(atoms) == 0:
            print("  [No atoms found for this species]\n")
            continue

        # select a few atoms to show
        show_atoms = atoms[:max_atoms]

        for atom in show_atoms:
            vec = atom_map[atom][:max_frames]   # show only first frames
            vec_str = "".join(str(x) for x in vec)
            print(f"  atom {atom:6d} : {vec_str}")

        print()  # blank line


# =====================================================================
# Main
# =====================================================================
if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage: python path_to_script oxy_table.csv")
        sys.exit(1)

    path = sys.argv[1]

    print("Loading and extracting species…")
    species_maps, n_frames = build_species_maps(path)

    print(f"Total frames detected: {n_frames}")
    for sp in TARGET_SPECIES:
        count_frames = len(species_maps[sp])
        print(f"  {sp}: seen in {count_frames} frames")

    print("\nBuilding indicator vectors…")
    indicators = build_indicator_vectors(species_maps, n_frames)

    preview_indicators(indicators)
