import pandas as pd
import numpy as np
from collections import defaultdict
from tqdm import tqdm


###############################################################
# --- 1. SPECIES MAPS (FINE / MIDDLE / COARSE) ---------------
###############################################################

# FINE: everything preserved
SPECIES_FINE = {
    "dissociated_OH" : "dissociated_OH",
    "remaining_OH"   : "remaining_OH",
    "vac_water_left" : "vac_water_left",
    "vac_water_right": "vac_water_right",
    "vac_water_up"   : "vac_water_up",
    "vac_water_down" : "vac_water_down",
    "OH_on_vacancy"  : "OH_on_vacancy",   # will be memory-resolved
}

# MIDDLE: unify left/right
SPECIES_MIDDLE = {
    "vac_water_left" : "vac_water_lr",
    "vac_water_right": "vac_water_lr",
    "vac_water_up"   : "vac_water_up",
    "vac_water_down" : "vac_water_down",
    "OH_on_vacancy"  : "vac_water_lr",    # after memory fix
}

# COARSE: unify all vacancy waters
SPECIES_COARSE = {
    "vac_water_left" : "vac_water",
    "vac_water_right": "vac_water",
    "vac_water_up"   : "vac_water",
    "vac_water_down" : "vac_water",
    "OH_on_vacancy"  : "vac_water",       # after memory fix
}


###############################################################
# --- 2. MAIN INDICATOR CONSTRUCTION ---------------------------
###############################################################

def build_all_indicator_vectors(csv_path):
    """
    Builds three separate indicator functions:
       fine   – fully resolved species
       mid    – unify left/right
       coarse – unify all vac waters

    Returns:
       indicators_fine, indicators_mid, indicators_coarse, n_frames
    """

    df = pd.read_csv(csv_path)

    # --- number of frames ---
    n_frames = df["frame"].max() + 1

    # --- persistent left/right memory for OH_on_vacancy ---
    persistent_orientation = {}   # atom idx → "left"/"right"

    # --- indicator dicts ---
    indicators_fine   = defaultdict(lambda: defaultdict(lambda: [0]*n_frames))
    indicators_mid    = defaultdict(lambda: defaultdict(lambda: [0]*n_frames))
    indicators_coarse = defaultdict(lambda: defaultdict(lambda: [0]*n_frames))

    # Iterate frame by frame (crucial for persistent orientation)
    for f in tqdm(range(n_frames), desc="building indicators"):

        df_f = df[df["frame"] == f]

        for _, row in df_f.iterrows():

            idx = int(row["atom_index"])
            sp  = row["species"]

            # ----------------------------------------------------------
            # 1) resolve OH_on_vacancy using persistent orientation
            # ----------------------------------------------------------
            if sp == "OH_on_vacancy":

                if idx in persistent_orientation:
                    side = persistent_orientation[idx]   # "left" / "right"
                else:
                    # infer from current row_tag (reliable)
                    rt = row["row_label"]
                    if rt == "vac_water_left" or rt == "left":
                        side = "left"
                    elif rt == "vac_water_right" or rt == "right":
                        side = "right"
                    else:
                        # fallback: guess from x value
                        side = "left" if row["atom_index"] % 2 == 0 else "right"

                # update memory
                persistent_orientation[idx] = side

                # transform species:
                sp = f"vac_water_{side}"

            # ----------------------------------------------------------
            # 2) Insert into fine indicator set
            # ----------------------------------------------------------
            sp_fine = SPECIES_FINE.get(sp)
            if sp_fine is not None:
                indicators_fine[sp_fine][idx][f] = 1

            # ----------------------------------------------------------
            # 3) Middle-grain indicators
            # ----------------------------------------------------------
            sp_mid = SPECIES_MIDDLE.get(sp)
            if sp_mid is not None:
                indicators_mid[sp_mid][idx][f] = 1

            # ----------------------------------------------------------
            # 4) Coarse indicators
            # ----------------------------------------------------------
            sp_coarse = SPECIES_COARSE.get(sp)
            if sp_coarse is not None:
                indicators_coarse[sp_coarse][idx][f] = 1

    return indicators_fine, indicators_mid, indicators_coarse, n_frames
