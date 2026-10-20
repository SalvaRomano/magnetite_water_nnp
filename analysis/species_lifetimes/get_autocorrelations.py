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
    "dissociated_OH" : "dissociated_OH",
    "remaining_OH"   : "remaining_OH",
    "vac_water_up"   : "vac_water_up",
    "vac_water_down" : "vac_water_down",
    "vac_water_left" : "vac_water",
    "vac_water_right": "vac_water",
    "OH_on_vacancy"  : "vac_water",
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
    return C 

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

    results_norm_then_avg = {}
    results_avg_then_norm = {}
    results_label = {}
  
    all_atom_indicators={}
    all_atom_acfs={}
    # =================================================================
    # Process each species (no CSV re-read)
    # =================================================================
    for sp in TARGET_SPECIES:
        print(f"\nProcessing species {sp}…")
        
        atom_indicators = {}
        atom_acfs = {}

           
        frame_dict = species_frames[sp]
        C_sum = np.zeros(tau_max + 1)
        C_norm_sum = np.zeros(tau_max + 1)
        n_eff = 0

        for atom, vec in tqdm(indicator_vectors_for_species(frame_dict, n_frames), desc=sp):
            if vec.sum() == 0:
                continue

            C = acf_indicator(vec, tau_max)
            #not normalized
            C_sum += C
            #normalized
            C_norm_sum += C/C[0]
            
            n_eff     += 1
            
            # Store indicator and raw ACF in dictionaries
            atom_indicators[atom] = vec.astype(np.uint8)   # compact
            atom_acfs[atom] = C



        if n_eff>0:
            C_norm_avg = C_norm_sum / n_eff
        else:
            C_norm_avg = np.zeros(tau_max+1)
        
        if C_sum[0] >0:
            C_avg_norm = C_sum/C_sum[0]
        else:
            C_avg_norm = np.zeros(tau_max+1)

        results_norm_then_avg[f"C_{sp}"] = C_norm_avg
        results_avg_then_norm[f"C_{sp}"] = C_avg_norm
        
        all_atom_indicators[sp]=atom_indicators
        all_atom_acfs[sp]=atom_acfs
        
       

    # =================================================================
    # Save results
    # =================================================================

    pd.DataFrame(results_norm_then_avg).to_csv("autocorr_norm_avg.csv", index=False)
    pd.DataFrame(results_avg_then_norm).to_csv("autocorr_avg_norm.csv", index=False)
    np.savez("atomic_indicators.npz", atomic_indicators=all_atom_indicators)
    np.savez("atomic_acfs.npz", atomic_acfs=all_atom_acfs)
   
    print("\nSaved autocorr.csv \natomic_indicators.npz \natomic_acfs.npz")

