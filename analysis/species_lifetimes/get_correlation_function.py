#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
from collections import defaultdict
from tqdm import tqdm
import tidynamics


# =====================================================================
# Species grouping rules
# =====================================================================

SPECIES_MAP = {
    "dissociated_OH": "dissociated_OH",
    "remaining_OH": "remaining_OH",
    "vac_water_up": "vac_water_up",
    "vac_water_down": "vac_water_down",
    "vac_water_left": "vac_water_side",
    "vac_water_right": "vac_water_side",
    "OH_on_vacancy"  : "OH_on_vacancy",
}

TARGET_SPECIES = sorted(set(SPECIES_MAP.values()))


# =====================================================================
# Step 1 — One pass: build frame → atoms structure
# =====================================================================
def build_species_frames(csv_path, chunksize=200000):
    """
    Returns:
        species_frames: {species: {frame: set(atom_ids)}}
        n_frames: total number of frames
    """
    species_frames = {sp: defaultdict(set) for sp in TARGET_SPECIES}
    max_frame = -1

    for chunk in pd.read_csv(
        csv_path, chunksize=chunksize,
        usecols=["frame", "atom_index", "species"]
    ):
        for _, row in chunk.iterrows():
            raw = row["species"]
            if raw not in SPECIES_MAP:
                continue

            sp = SPECIES_MAP[raw]
            f = int(row["frame"])
            a = int(row["atom_index"])
            species_frames[sp][f].add(a)
            max_frame = max(max_frame, f)

    n_frames = max_frame + 1
    return species_frames, n_frames


# =====================================================================
# Build indicator vectors from species_frames
# =====================================================================
def indicator_vectors_for_species(frame_dict, n_frames):
    """
    Yields (atom, indicator_vector) for the given species.
    Only atoms that ever appear are included.
    """
    if not frame_dict:
        return

    all_atoms = sorted(set().union(*frame_dict.values()))

    for atom in all_atoms:
        vec = np.zeros(n_frames, dtype=bool)
        for f, atoms in frame_dict.items():
            if atom in atoms:
                vec[f] = True
        yield atom, vec


# =====================================================================
# ACF using tidynamics
# =====================================================================
def acf_indicator(vec, tau_max=None):
    vec = np.asarray(vec, dtype=float)
    C = tidynamics.acf(vec)
    if tau_max is not None:
        C = C[:tau_max + 1]
    return C / C[0]


# =====================================================================
# MAIN
# =====================================================================
if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Fast one-pass ACF calculation for surface OH species")
    parser.add_argument("oxy_table")
    parser.add_argument("--tau-max", type=int, default=None)
    args = parser.parse_args()

    csv_path = args.oxy_table

    print("\nReading CSV in one pass to build species_frames…")
    species_frames, n_frames = build_species_frames(csv_path)
    tau_max = args.tau_max or (n_frames - 1)
    print(f"  Found {n_frames} frames")
    print(f"  tau_max = {tau_max}")

    results = {}
    lifetimes = {}

    # =================================================================
    # Process each species (no CSV re-read)
    # =================================================================
    for sp in TARGET_SPECIES:
        print(f"\nProcessing species {sp}…")

        frame_dict = species_frames[sp]
        C_sum = np.zeros(tau_max + 1)
        n_eff = 0

        for atom, vec in tqdm(indicator_vectors_for_species(frame_dict, n_frames), desc=sp):
            if vec.sum() == 0:
                continue

            C = acf_indicator(vec, tau_max)
            C_sum += C
            n_eff += 1

        if n_eff > 0:
            C_avg = C_sum / n_eff
        else:
            C_avg = np.zeros(tau_max + 1)

        results[f"C_{sp}"] = C_avg

        # Lifetime
        dt_ps = 0.004
        tau_ps = float(np.sum(C_avg) * dt_ps)
        lifetimes[sp] = tau_ps
        print(f"  Lifetime({sp}) = {tau_ps:.6f} ps")


    # =================================================================
    # Save results
    # =================================================================

    pd.DataFrame(results).to_csv("autocorr.csv", index=False)
    print("\nSaved autocorr.csv")

    pd.DataFrame({
        "species": list(lifetimes.keys()),
        "lifetime_ps": list(lifetimes.values())
    }).to_csv("lifetimes.csv", index=False)

    print("Saved lifetimes.csv\nDone.")
