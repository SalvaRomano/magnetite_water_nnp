#!/usr/bin/env python3

import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import defaultdict
import tidynamics

from build_indicators import build_all_indicator_vectors


# ---------------------------------------------------------------------
# ACF wrapper
# ---------------------------------------------------------------------
def acf_indicator(vec, tau_max=None):
    vec = np.asarray(vec, dtype=float)
    C = tidynamics.acf(vec)   # no mean subtraction (correct for indicators)
    if tau_max is not None:
        C = C[:tau_max + 1]
    return C


# ---------------------------------------------------------------------
# Process indicator dict: compute ACFs, averaged ACFs, normalized ACFs
# ---------------------------------------------------------------------
def process_indicator_set(indicators, tau_max):

    results_normavg = {}   # <C_i/C_i(0)>_atoms
    results_avgnorm = {}   # C_total / C_total(0)

    all_indicators = {}
    all_acfs = {}

    for sp, atom_vectors in indicators.items():
        print(f"\nProcessing species: {sp}")

        C_sum = np.zeros(tau_max + 1)
        C_norm_sum = np.zeros(tau_max + 1)
        n_eff = 0

        all_indicators[sp] = {}
        all_acfs[sp] = {}

        for atom, vec in tqdm(atom_vectors.items(), desc=sp):
            vec = np.asarray(vec, dtype=float)
            if vec.sum() == 0:
                continue

            C = acf_indicator(vec, tau_max)

            # store raw
            all_indicators[sp][atom] = vec.astype(np.uint8)
            all_acfs[sp][atom] = C

            # accumulations
            C_sum      += C
            C_norm_sum += C / C[0]
            n_eff      += 1

        # averages
        if n_eff > 0:
            C_normavg = C_norm_sum / n_eff     # < C_i(τ)/C_i(0) >
        else:
            C_normavg = np.zeros(tau_max + 1)

        if C_sum[0] > 0:
            C_avgnorm = C_sum / C_sum[0]       # C_total(τ)/C_total(0)
        else:
            C_avgnorm = np.zeros(tau_max + 1)

        results_normavg[f"C_{sp}"] = C_normavg
        results_avgnorm[f"C_{sp}"] = C_avgnorm

    return results_normavg, results_avgnorm, all_indicators, all_acfs


# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------
if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("oxy_table")
    parser.add_argument("--tau-max", type=int, default=None)
    args = parser.parse_args()

    csv_path = args.oxy_table

    print("\n=== Building fine, middle, and coarse indicator sets ===")
    indicators_fine, indicators_mid, indicators_coarse, n_frames = \
        build_all_indicator_vectors(csv_path)

    tau_max = args.tau_max or (n_frames - 1)
    print(f"  Found {n_frames} frames")
    print(f"  tau_max = {tau_max}\n")

    # master containers (merged across all resolutions)
    merged_normavg = {}
    merged_avgnorm = {}
    merged_indicators = {}
    merged_acfs = {}

    # ------------------------------------------------------------
    # Run all three resolutions in order: FINE, MID, COARSE
    # ------------------------------------------------------------
    for label, inds in [
        ("fine",   indicators_fine),
        ("middle", indicators_mid),
        ("coarse", indicators_coarse)
    ]:

        print(f"\n=== Processing {label} indicator set ===")

        normavg, avgnorm, ind, acf = process_indicator_set(inds, tau_max)

        # merge results
        merged_normavg.update(normavg)
        merged_avgnorm.update(avgnorm)
        merged_indicators.update(ind)
        merged_acfs.update(acf)

    # ------------------------------------------------------------
    # Save final merged files
    # ------------------------------------------------------------

    pd.DataFrame(merged_normavg).to_csv("all_autocorr_normavg.csv", index=False)
    pd.DataFrame(merged_avgnorm).to_csv("all_autocorr_avgnorm.csv", index=False)

    np.savez("all_atomic_indicators.npz", atomic_indicators=merged_indicators)
    np.savez("all_atomic_acfs.npz",        atomic_acfs=merged_acfs)

    print("\n=== Saved final merged output ===")
    print("  all_autocorr_normavg.csv")
    print("  all_autocorr_avgnorm.csv")
    print("  all_atomic_indicators.npz")
    print("  all_atomic_acfs.npz\n")
