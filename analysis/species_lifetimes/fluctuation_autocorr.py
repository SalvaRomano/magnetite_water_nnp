#!/usr/bin/env python3

import numpy as np
import pandas as pd


# -----------------------------
# 1. Load atomic ACFs
# -----------------------------
# atomic_acfs.npz was saved as:
# np.savez("atomic_acfs.npz", atomic_acfs=all_atomic_acfs)
# where all_atomic_acfs[species][atom] = C_raw(τ)

acf_data = np.load("atomic_acfs.npz", allow_pickle=True)
all_atomic_acfs = acf_data["atomic_acfs"].item()  # nested dict: sp -> {atom: acf_array}


# -----------------------------
# 2. Compute p_M1 and p_M2 per species
# -----------------------------
p_M1 = {}  # simple mean of p_i
p_M2 = {}  # Σ p_i^2 / Σ p_i

for sp, acf_dict in all_atomic_acfs.items():
    # p_i = C_i^raw(0)
    p_i = np.array([acf[0] for acf in acf_dict.values()], dtype=float)

    if len(p_i) == 0 or np.all(p_i == 0):
        p_M1[sp] = 0.0
        p_M2[sp] = 0.0
        continue

    # Method 1 plateau: average p_i
    p_M1[sp] = float(p_i.mean())

    # Method 2 plateau: Σ p_i^2 / Σ p_i
    p_M2[sp] = float( (p_i**2).sum() / p_i.sum() )


# Save plateaus to a small table
plateau_df = pd.DataFrame({
    "species": list(p_M1.keys()),
    "p_M1": [p_M1[sp] for sp in p_M1],
    "p_M2": [p_M2[sp] for sp in p_M2],
})
plateau_df.to_csv("plateaus_M1_M2.csv", index=False)


# -----------------------------
# 3. Load averaged ACFs (M1, M2)
# -----------------------------
df_M1 = pd.read_csv("autocorr_norm_avg.csv")
df_M2 = pd.read_csv("autocorr_avg_norm.csv")

# We expect column names like "C_M1_<species>" and "C_M2_<species>"
cols_M1 = [c for c in df_M1.columns if c.startswith("C_")]
cols_M2 = [c for c in df_M2.columns if c.startswith("C_")]

# -----------------------------
# 4. Build fluctuation ACFs: (C - p)/(1 - p)
# -----------------------------
df_M1_fluct = df_M1.copy()
df_M2_fluct = df_M2.copy()

for col in cols_M1:
    sp = col.replace("C_M1_", "")
    p = p_M1.get(sp, 0.0)
    if p < 1.0:
        df_M1_fluct[col] = (df_M1[col] - p) / (1.0 - p)
    else:
        # pathological case: p=1, everything always in state
        df_M1_fluct[col] = 0.0

for col in cols_M2:
    sp = col.replace("C_M2_", "")
    p = p_M2.get(sp, 0.0)
    if p < 1.0:
        df_M2_fluct[col] = (df_M2[col] - p) / (1.0 - p)
    else:
        df_M2_fluct[col] = 0.0


# -----------------------------
# 5. Save the fluctuation ACFs
# -----------------------------
df_M1_fluct.to_csv("autocorr_M1_fluct.csv", index=False)
df_M2_fluct.to_csv("autocorr_M2_fluct.csv", index=False)

print("Saved:")
print("  plateaus_M1_M2.csv")
print("  autocorr_M1_fluct.csv")
print("  autocorr_M2_fluct.csv")
