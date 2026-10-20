#!/usr/bin/env python3

import numpy as np
from math import sqrt
from tqdm import tqdm
import MDAnalysis as mda
import argparse

from scipy.spatial import cKDTree


def read_first_frame(dumpfile):
    """
    Read only the first frame of a LAMMPS dump to get atom types.
    Returns:
        n_atoms (int),
        type_nums (np.ndarray of ints)
    """
    atom_lines = []
    with open(dumpfile, "r", errors="ignore") as f:
        for line in f:
            if line.startswith("ITEM: ATOMS"):
                header = line.strip().split()[2:]
                if "type" not in header:
                    raise ValueError("LAMMPS dump does not contain a 'type' column.")
                type_col = header.index("type")
                for l in f:
                    if l.startswith("ITEM:"):
                        break
                    atom_lines.append(l.strip().split()[type_col])
                break

    # convert to int numpy array
    type_nums = np.array(atom_lines, dtype=int)
    n_atoms = len(type_nums)
    return n_atoms, type_nums

def build_topo(dumpfile, atom_names):
    # --- Build minimal topology directly from dump ---
    n_atoms, type_nums = read_first_frame(dumpfile)
    
    type2name = {1: atom_names[0], 2: atom_names[1], 3: atom_names[2]}
    names = np.array([type2name.get(t, f"T{t}") for t in type_nums])
    ids = np.arange(1, n_atoms + 1)
    
    u = mda.Universe.empty(n_atoms, trajectory=True)
    u.add_TopologyAttr("types", type_nums)   # <── added line
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("ids", ids)
    u.load_new(dumpfile, format="LAMMPSDUMP")
    return u
    
  




# ---------- clustering utilities ----------
# ---------- Fe layer identification ----------
def locate_fe_layers(fe_atoms, fe_buffer,
                     tet_offset=1.0, tet_z_tol=0.6):
    """
    Identify octahedral surface Fe and tetrahedral subsurface Fe layers.

    Inputs
    ------
    fe_atoms : AtomGroup or array-like of Fe atoms (has positions and indices)
    fe_buffer : float
        thickness used to capture the outermost octahedral Fe
    tet_offset : float
        expected vertical distance between octahedral and tetrahedral Fe layers
    tet_z_tol : float
        selection tolerance for identifying tetrahedral Fe layer
    
    Returns
    -------
    fe_top_sel : Fe atoms in the top octahedral layer
    fe_bot_sel : Fe atoms in the bottom octahedral layer
    z_top : float, mean z of top octahedral layer
    z_bot : float, mean z of bottom octahedral layer
    
    fe_tet_top_sel : Fe atoms in the subsurface tetrahedral layer (top)
    fe_tet_bot_sel : Fe atoms in the subsurface tetrahedral layer (bottom)
    z_tet_top : float, mean z of top tetrahedral layer
    z_tet_bot : float, mean z of bottom tetrahedral layer
    """

    zvals = fe_atoms.positions[:, 2]
    zmax, zmin = zvals.max(), zvals.min()

    # ---- Surface (octahedral) Fe selection ----
    fe_top_sel = fe_atoms[zvals >= zmax - fe_buffer]
    fe_bot_sel = fe_atoms[zvals <= zmin + fe_buffer]

    # compute surface z-centers
    z_top = float(np.mean(fe_top_sel.positions[:, 2])) if len(fe_top_sel) > 0 else np.nan
    z_bot = float(np.mean(fe_bot_sel.positions[:, 2])) if len(fe_bot_sel) > 0 else np.nan

    # ---- Subsurface tetrahedral Fe selection ----
    z_tet_top_center = z_top - tet_offset   # 1 Å below top surface
    z_tet_bot_center = z_bot + tet_offset   # 1 Å above bottom surface

    mask_tet_top = (
        (zvals > z_tet_top_center - tet_z_tol) &
        (zvals < z_tet_top_center + tet_z_tol)
    )
    mask_tet_bot = (
        (zvals > z_tet_bot_center - tet_z_tol) &
        (zvals < z_tet_bot_center + tet_z_tol)
    )

    fe_tet_top_sel = fe_atoms[mask_tet_top]
    fe_tet_bot_sel = fe_atoms[mask_tet_bot]

    # compute tet-layer mean z values
    z_tet_top = float(np.mean(fe_tet_top_sel.positions[:, 2])) if len(fe_tet_top_sel) > 0 else np.nan
    z_tet_bot = float(np.mean(fe_tet_bot_sel.positions[:, 2])) if len(fe_tet_bot_sel) > 0 else np.nan

    return (fe_top_sel, fe_bot_sel, z_top, z_bot,
            fe_tet_top_sel, fe_tet_bot_sel, z_tet_top, z_tet_bot)

    

def projections(atoms, normal):
    return atoms.positions[:, 0]*normal[0]+atoms.positions[:, 1]*normal[1]+atoms.positions[:, 2]*normal[2]
        

def cluster_1d(values, row_cluster_threshold):
    """Simple 1D clustering by grouping sorted values with small gaps."""
    if len(values) == 0:
        return []
    order = np.argsort(values)
    vals = np.asarray(values)[order]
    clusters, current = [], [order[0]]
    
    for i in range(1, len(vals)):
        if (vals[i] - vals[i - 1]) <= row_cluster_threshold:
            current.append(order[i])
        else:
            clusters.append(current)
            current = [order[i]]
    clusters.append(current)
    return [(float(np.mean(np.asarray(values)[inds])), inds) for inds in clusters]


def pbc_distance(coord1, coord2, box):
    cell_lengths=box[:3]
    """
    Calculates the distance between two coordinates with periodic boundary conditions.

    Args:
        coord1 (np.array): The first coordinate (e.g., [x1, y1, z1]).
        coord2 (np.array): The second coordinate (e.g., [x2, y2, z2]).
        cell_lengths (np.array): The lengths of the simulation box in each dimension (e.g., [Lx, Ly, Lz]).

    Returns:
        float: The periodic distance.
    """
    delta = coord1 - coord2
    # Apply minimum image convention
    for i in range(len(cell_lengths)):
        if abs(delta[i]) > cell_lengths[i] * 0.5:
            delta[i] = delta[i] - cell_lengths[i] * np.sign(delta[i])
    
    return np.sqrt(np.sum(delta**2))

def create_PBC_image(h_pos, h_idx, box):
    Lx, Ly, Lz = box[:3]
    # --- build periodic hydrogen images (only in xy plane for a slab) ---
    
    shifts = np.array([
        [0, 0, 0],
        [Lx, 0, 0], [-Lx, 0, 0],
        [0, Ly, 0], [0, -Ly, 0],
        [Lx, Ly, 0], [Lx, -Ly, 0],
        [-Lx, Ly, 0], [-Lx, -Ly, 0]
    ])
    
    h_images = []
    h_indices = []
    for s in shifts:
        h_images.append(h_pos + s)
        h_indices.append(h_idx)     # same IDs for each image
    h_all = np.vstack(h_images)
    h_idx_all = np.concatenate(h_indices)
    return h_idx_all, h_all



def row_label_for_projection(proj_val, centers, row_spacing,  onrow_fraction, L_diag):
    
    """Return 'oct' if within 1/8 spacing (periodic along diagonal)."""

    
    if centers is None or len(centers) == 0:
        return "unknown"
    # diagonal-periodic distance
    d = np.abs(proj_val - centers)
    d = np.minimum(d, L_diag - d)  # PBC across diagonal
    if np.min(d) < onrow_fraction * row_spacing:
        return "oct"
    return "tet"




#oxygen classification
def classify_surface_oxygens(h_atoms, o_atoms , z_top, z_bot,
                             surf_bounds, layer_bounds,
                             frame, file_top, file_bot, f_sum_top, f_sum_bot, dipole_theta_boundary,
                             h_cutoff=1.2,
                             fe_top_centers=None,
                             fe_bot_centers=None,
                             box=None,
                             onrow_fraction=0.125,
                             vac_cutoff=3,
                             wat_on_vac_h_cutoff=1.3,
                             vac_top=None,
                             vac_bot=None):
    """
    Identify O atoms in surface and first-layer regions for both surfaces.
    For each O:
      – counts nearby H atoms (KDTree, within h_cutoff)
      – adds label 'oct' or 'tet' using Fe-row centers
      
    Columns:
        frame,atom_index,surface,region,h_count,h_indices,row_label
    """
    # Lateral box lengths for diagonal-periodic projection
    Lx, Ly, _ = box[:3]
    L_diag = np.sqrt(Lx**2 + Ly**2)

    z_o_all = o_atoms.positions[:, 2]
    z_h_all = h_atoms.positions[:, 2]
    idx_o_all = o_atoms.indices
    idx_h_all = h_atoms.indices

    # --- region masks for O ---
    o_mask_surface_top = ((z_o_all < z_top + surf_bounds[1]) &
                          (z_o_all > z_top + surf_bounds[0]))
    o_mask_surface_bot = ((z_o_all > z_bot - surf_bounds[1]) &
                          (z_o_all < z_bot - surf_bounds[0]))
    o_mask_layer_top = ((z_o_all < z_top + layer_bounds[1]) &
                        (z_o_all >= z_top + layer_bounds[0]))
    o_mask_layer_bot = ((z_o_all > z_bot - layer_bounds[1]) &
                        (z_o_all <= z_bot - layer_bounds[0]))
    
    # --- hydrogens near surfaces ---
    h_mask_up = ((z_h_all < z_top + layer_bounds[1] + h_cutoff) &
                 (z_h_all > z_top - h_cutoff))
    h_mask_down = ((z_h_all > z_bot - layer_bounds[1] - h_cutoff) &
                   (z_h_all < z_bot + h_cutoff))

    h_up_pos, h_up_idx = h_atoms.positions[h_mask_up], idx_h_all[h_mask_up]
    h_down_pos, h_down_idx = h_atoms.positions[h_mask_down], idx_h_all[h_mask_down]

    

    # --- process all four regions ---
    top_normal = [1/np.sqrt(2), 1/np.sqrt(2), 0]
    bottom_normal = [1/np.sqrt(2), -1/np.sqrt(2), 0]

    
   
    
    # --- common KDTree writer with row label ---
    def write_region(file_handle, o_mask, surface, region, h_pos, h_idx,
                 fe_centers, normal, box, vac_cutoff, onrow_fraction, vac_positions,  wat_on_vac_h_cutoff , dipole_theta_boundary):
        """
        Write one surface/region block.
        For 'layer' regions, label the vac_cutoff tet oxygens as 'vacancy_water'.
        """
       

        o_pos = o_atoms.positions[o_mask]
        o_idx = idx_o_all[o_mask]
        o_z   = z_o_all[o_mask]
        if len(o_pos) == 0:
            return

        # --- PBC hydrogen images ---
        h_idx_all, h_all = create_PBC_image(h_pos, h_idx, box)

        # --- KDTree search ---
        tree = cKDTree(h_all)
        neighbours = tree.query_ball_point(o_pos, h_cutoff)

        # --- Fe-row geometry ---
        if fe_centers is not None and len(fe_centers) > 1:
            row_spacing = np.mean(np.diff(np.sort(fe_centers)))
        else:
            row_spacing = np.inf

        o_proj = projections(o_atoms[o_mask], normal)

        # --- build temporary records --- 
        records = []
        for oi, nlist, proj, zval, pos in zip(o_idx, neighbours, o_proj, o_z, o_pos):
            if len(nlist) == 0:
                hcount, hlist_str = 0, "[]"
            else:
                h_real = [int(h_idx_all[i]) for i in nlist]
                hcount = len(h_real)
                hlist_str = "[" + ", ".join(map(str, h_real)) + "]"

            row_tag = row_label_for_projection(proj, fe_centers, row_spacing,
                                               onrow_fraction, L_diag)
            records.append(dict(idx=oi,
                                hcount=hcount,
                                hlist_str=hlist_str,
                                row_tag=row_tag,
                                z=zval,
                                pos=pos))

    
        # --- vacancy-based refinement for tet waters in first layer ---
        if region == "layer" and vac_positions is not None and len(vac_positions) > 0:
        
            # indices of tet-water records
            tet_indices = [i for i, r in enumerate(records) if r["row_tag"] == "tet"]
            if not tet_indices:
                pass
            else:
                # Extract all tet positions
                tet_pos = np.array([records[i]["pos"] for i in tet_indices])
                tet_idx_arr = np.array(tet_indices)

                # ---- Step 1: select only tet oxygens near borders ----
                R = vac_cutoff      # vacancy water cutoff
                x = tet_pos[:,0]
                y = tet_pos[:,1]

                near_x = (x < R) | (x > Lx - R)
                near_y = (y < R) | (y > Ly - R)
                near_edge = near_x | near_y

                tet_edge_pos = tet_pos[near_edge]
                tet_edge_idx = tet_idx_arr[near_edge]

                # ---- Step 2: PBC replicate only tet on border using existing util ----
                tet_edge_idx_all, tet_edge_pos_all = create_PBC_image(
                    tet_edge_pos, tet_edge_idx, box
                )

                # interior (non-replicated) tet positions
                tet_interior_pos = tet_pos[~near_edge]
                tet_interior_idx = tet_idx_arr[~near_edge]

                # ---- build full KDTree set ----
                tet_all_pos = np.vstack([tet_interior_pos, tet_edge_pos_all])
                tet_all_idx = np.concatenate([tet_interior_idx, tet_edge_idx_all])

                # ---- build KDTree ----
                tree_tet = cKDTree(tet_all_pos)

                used = set()

                # ---- Step 3: loop over vacancies ----
                for v in vac_positions:   # each vac is a 3D (x,y,z)
                    # list of neighbor tet oxygens within R_vac
                    nbrs = tree_tet.query_ball_point(v, r=R)
                    if not nbrs:
                        continue
                    
                    # choose closest tet oxygen
                    best = None
                    best_d2 = None
                    for j in nbrs:
                        oi = tet_all_idx[j]   # original record index
                        if oi in used:
                            continue
                        
                        p = tet_all_pos[j]
                        dx = p[0] - v[0]
                        dy = p[1] - v[1]
                        dz = p[2] - v[2]
                        d2 = dx*dx + dy*dy + dz*dz

                        if best_d2 is None or d2 < best_d2:
                            best_d2 = d2
                            best = oi

                    if best is not None:
                        records[best]["row_tag"] = "on_vacancy"
                        used.add(best)

        #dipole analysis
   
        for r in records:
            r["orientation_tag"] = ""         
            # ----- Unified vacancy refinement: H-count + dipole + orientation -----
            if r["row_tag"] == "on_vacancy":
            
                # --- KDTree query: find hydrogens near this oxygen with greater cutoff---
                neighbors = tree.query_ball_point(r["pos"], wat_on_vac_h_cutoff)

                # collect OH vectors relative to O (dx,dy,dz)
                Hs = []
                for idxh in neighbors:
                    hpos = h_all[idxh]
                    # displacement vector: hpos is already PBC-corrected image
                    dx, dy, dz = hpos - r["pos"]
                    Hs.append(np.array([dx, dy, dz]))

                # ---- vacancy-tailored hcount ----
                r["hcount"] = len(Hs)
                if len(Hs) == 2:                
                    # ---- dipole from two hydrogens ----
                    v1, v2 = Hs
                    dip = v1 + v2

                    dn = np.linalg.norm(dip)
                   
                    dunit = dip / dn
                   
                    mux, muy, muz = dunit

                    # ---- choose surface normal ----
                    if surface == "top":
                        nvec = np.array([0.0, 0.0, 1.0])
                    else:
                        nvec = np.array([0.0, 0.0, -1.0])

                    # ---- angles ----
                    dot = float(np.dot(dunit, nvec))
                    theta = np.degrees(np.arccos(dot))
                    phi   = np.degrees(np.arctan2(muy, mux))

                    # ---- orientation class ----
                    if theta < dipole_theta_boundary[0]:
                        orient = "_up"
                    elif theta < dipole_theta_boundary[1]:
                        orient = "_left" if phi > 0 else "_right"
                    else:
                        orient = "_down"

                    # store
                    r["orientation_tag"] = orient


        # --- counters for summary ---
        counts = {
            "dissociated_OH": 0,
            "remaining_OH": 0,
            "oct_water": 0,
            "tet_water": 0,
            "vac_water": 0,
            "OH_on_vacancy": 0,
            "hydronium_on_vacancy": 0,
            "degenerate": 0,
            "vac_up":0,
            "vac_down":0,
            "vac_middle":0,
            "0_diss": 0,
            "1_diss": 0,
            "2_diss": 0,
        }
        
        
        # ---------------------------------------------------------
        # --- vacancy-based dissociation counting via OH oxygens ---
        # ---------------------------------------------------------
        vac_H0 = vac_H1 = vac_H2 = 0
        
        if region == "surface":
        
            # collect dissociated OH positions
            diss_O_pos = []
            for r in records:
                if r["hcount"] == 1:
                    diss_O_pos.append(r["pos"])
        
            if len(diss_O_pos) == 0:
                vac_H0 = len(vac_positions)
            else:
                # --- PBC-expansion for OH oxygens ---
                diss_O_pos = np.array(diss_O_pos)
                dummy_idx = np.arange(len(diss_O_pos))
                diss_idx_all, diss_pos_all = create_PBC_image(diss_O_pos, dummy_idx, box)
        
                # KDTree on expanded OH oxygens
                tree_diss = cKDTree(diss_pos_all)
        
                R_OOH = 3.5  # Å cutoff
        
                for v in vac_positions:
                    nbrs = tree_diss.query_ball_point(v, r=R_OOH)
                    k = len(nbrs)
        
                    if k == 0:
                        vac_H0 += 1
                    elif k == 1:
                        vac_H1 += 1
                    elif k == 2:
                        vac_H2 += 1
        
                    

        counts["0_diss"] = vac_H0
        counts["1_diss"] = vac_H1
        counts["2_diss"] = vac_H2
        
        # --- write all records ---
        for r in records:
            # skip surface oxygens with no hydrogens attached
            if region == "surface" and r["hcount"] == 0:
                continue
                        
            # --- determine species label ---
            species = "degenerate"  # default fallback

            if region == "surface" and r["hcount"] == 1:
                species = "dissociated_OH"

            elif region == "layer" and r["hcount"] == 1:
                if r["row_tag"] == "oct":
                    species = "remaining_OH"
                elif r["row_tag"] == "on_vacancy":
                    # rare: OH sitting in vacancy pocket
                    species = "OH_on_vacancy"

            elif region == "layer" and r["hcount"] == 2:
                if r["row_tag"] == "oct":
                    species = "oct_water"
                elif r["row_tag"] == "tet":
                    species = "tet_water"
                elif r["row_tag"] == "on_vacancy":
                    species = "vac_water"

            elif region == "layer" and r["hcount"] == 3 and r["row_tag"] == "on_vacancy":
                    species = "hydronium_on_vacancy"
                    
                    
            #call back tet row tag
            if r["row_tag"] == "on_vacancy":
                tag = r["orientation_tag"]
                
                if tag == "_up":
                    counts["vac_up"] += 1

                elif tag == "_down":
                    counts["vac_down"] += 1
                elif ( tag == "_right" or tag == "_left" ): 
                    counts["vac_middle"] += 1
                #give again row_tag the corret meaning: "tet"   
                r["row_tag"] = "tet"
            
            #_update count
            counts[species] = counts.get(species, 0) + 1
            

            # --- write to detailed file ---
            file_handle.write(
                f'{frame},{r["idx"]},{region},{r["hcount"]},"{r["hlist_str"]}",{r["row_tag"]},{species}{r["orientation_tag"]}\n'
            )
        # --- after writing all records, return counts for this region ---
        return counts
    
    # gather counts from both regions for each surface
    counts_top = {"dissociated_OH":0,"0_diss":0,"1_diss":0,"2_diss":0,"remaining_OH":0,"vac_water":0, "vac_up":0, "vac_middle":0,
    "vac_down":0, "oct_water":0,"tet_water":0,"OH_on_vacancy":0,"hydronium_on_vacancy":0,"degenerate":0}
    counts_bot = {k:0 for k in counts_top}

    # ---- TOP ----
    c_top_surface = write_region(file_top, o_mask_surface_top, "top", "surface", h_up_pos, h_up_idx,
                                 fe_top_centers, top_normal, box, vac_cutoff, onrow_fraction, 
                                 vac_positions=vac_top, wat_on_vac_h_cutoff=wat_on_vac_h_cutoff, dipole_theta_boundary=dipole_theta_boundary )
    c_top_layer   = write_region(file_top,o_mask_layer_top, "top", "layer", h_up_pos, h_up_idx,
                                 fe_top_centers, top_normal, box, vac_cutoff, onrow_fraction, 
                                 vac_positions=vac_top, wat_on_vac_h_cutoff=wat_on_vac_h_cutoff, dipole_theta_boundary=dipole_theta_boundary)

    # ---- BOTTOM ----
    c_bot_surface = write_region(file_bot,o_mask_surface_bot, "bot", "surface", h_down_pos, h_down_idx,
                                 fe_bot_centers, bottom_normal, box, vac_cutoff, onrow_fraction, 
                                 vac_positions=vac_bot, wat_on_vac_h_cutoff=wat_on_vac_h_cutoff, dipole_theta_boundary=dipole_theta_boundary)
    c_bot_layer   = write_region(file_bot,o_mask_layer_bot, "bot", "layer", h_down_pos, h_down_idx,
                                 fe_bot_centers, bottom_normal, box, vac_cutoff, onrow_fraction,
                                 vac_positions=vac_bot, wat_on_vac_h_cutoff=wat_on_vac_h_cutoff, dipole_theta_boundary=dipole_theta_boundary)

    # combine counts
    for k in counts_top:
        counts_top[k] = c_top_surface.get(k,0) + c_top_layer.get(k,0)
        counts_bot[k] = c_bot_surface.get(k,0) + c_bot_layer.get(k,0)

    # --- write one summary line per surface ---
  
    f_sum_top.write(
        f'{frame},'
        f'{counts_top["dissociated_OH"]},{counts_top["0_diss"]},{counts_top["1_diss"]},{counts_top["2_diss"]},{counts_top["remaining_OH"]},{counts_top["vac_water"]},{counts_top["vac_up"]},{counts_top["vac_middle"]},{counts_top["vac_down"]},'
        f'{counts_top["oct_water"]},{counts_top["tet_water"]},{counts_top["OH_on_vacancy"]},{counts_top["hydronium_on_vacancy"]},{counts_top["degenerate"]}\n'
    )
    f_sum_bot.write(
        f'{frame},'
        f'{counts_bot["dissociated_OH"]},{counts_bot["0_diss"]},{counts_bot["1_diss"]},{counts_bot["2_diss"]},{counts_bot["remaining_OH"]},{counts_bot["vac_water"]},{counts_bot["vac_up"]},{counts_bot["vac_middle"]},{counts_bot["vac_down"]},'
        f'{counts_bot["oct_water"]},{counts_bot["tet_water"]},{counts_bot["OH_on_vacancy"]},{counts_bot["hydronium_on_vacancy"]},{counts_bot["degenerate"]}\n'
    )




def analyze_fe_rows_for_surface(
        fe_oct_sel,
        fe_tet_sel,
        row_cluster_threshold,
        row_normal,
        box
    ):
    """
    Analyze Fe rows for ONE surface (top or bottom), including:
      - octahedral surface rows
      - tetrahedral subsurface rows
      - PBC-aware merging of tet clusters into physical rows (for vacancy detection)
      - vacancy detection using correct tet row centers (not artificial merged centers)

    Returns:
        z_oct, fe_oct_centers,
        z_tet, fe_tet_centers,   # fe_tet_centers = true tet row centers (cluster_1d)
        vacancies  (array Nx3)
    """

    Lx, Ly, Lz = box[:3]

    positions_oct = fe_oct_sel.positions
    positions_tet = fe_tet_sel.positions

    # --- z-locations ---
    z_oct = float(np.mean(positions_oct[:, 2])) if len(positions_oct) else np.nan
    z_tet = float(np.mean(positions_tet[:, 2])) if len(positions_tet) else np.nan

    # --- normalize row normal ---
    row_normal = np.asarray(row_normal, float)
    row_normal /= np.linalg.norm(row_normal)

    nx, ny = row_normal[0], row_normal[1]

    # --- row direction in xy-plane (±45°) ---
    # Convention:
    #   ny >= 0  → anti-diagonal  (1, -1)
    #   ny <  0  → diagonal       (1,  1)
    if ny >= 0:
        row_dir_xy = np.array([1.0, -1.0])
    else:
        row_dir_xy = np.array([1.0,  1.0])
    row_dir_xy /= np.linalg.norm(row_dir_xy)
    row_dir = np.array([row_dir_xy[0], row_dir_xy[1], 0.0])

    # --- diagonal period ---
    L_diag = np.sqrt(Lx*Lx + Ly*Ly)
    half_diag = 0.5 * L_diag
    merge_tol = 0.25 * L_diag

    # -------------------------
    # 1) OCTAHEDRAL ROWS
    # -------------------------
    if len(positions_oct) > 0:
        oct_proj = positions_oct[:,0]*row_normal[0] + positions_oct[:,1]*row_normal[1]
        oct_clusters = cluster_1d(oct_proj, row_cluster_threshold)
        fe_oct_centers = [float(c) for c,_ in sorted(oct_clusters, key=lambda x: x[0])]
    else:
        fe_oct_centers = []

    # -------------------------
    # 2) TETRAHEDRAL ROWS
    # -------------------------
    vacancies = []

    if len(positions_tet) > 0:
        tet_proj = positions_tet[:,0]*row_normal[0] + positions_tet[:,1]*row_normal[1]
        tet_clusters = cluster_1d(tet_proj, row_cluster_threshold)
        if len(tet_clusters) == 0:
            fe_tet_centers = []
            return z_oct, fe_oct_centers, z_tet, fe_tet_centers, np.zeros((0,3))

        # sort clusters by their transverse center
        tet_clusters_sorted = sorted(tet_clusters, key=lambda x: x[0])
        centers = np.array([c for c,_ in tet_clusters_sorted])      # TRUE tet row centers
        fe_tet_centers = sorted(centers.tolist())                    # output these

        # for merging (to detect vacancies across PBC), work modulo L_diag
        centers_mod = centers % L_diag
        ncl = len(centers_mod)

        used = np.zeros(ncl, dtype=bool)
        physical_rows = []  # each entry is list of cluster indices that form one physical row

        # --- build physical row groups ---
        for i in range(ncl):
            if used[i]:
                continue
            ci = centers_mod[i]
            best_j = None
            best_err = 1e9
            for j in range(ncl):
                if i == j or used[j]:
                    continue
                cj = centers_mod[j]
                d = (cj - ci) % L_diag
                err = abs(d - half_diag)
                if err < best_err:
                    best_err = err
                    best_j = j

            if best_j is not None and best_err < merge_tol:
                physical_rows.append([i, best_j])
                used[i] = True
                used[best_j] = True
            else:
                physical_rows.append([i])
                used[i] = True
     
        # -------------------------
        # 3) PROCESS EACH PHYSICAL ROW
        # -------------------------
        vacancies = []
        for group in physical_rows:
         
            # collect segments for this physical row
            row_long_list = []
            row_pos_list  = []
         
           
            # for each cluster index in this physical row
            # we decide later which one to shift by +half_diag
            # but for now we build their base longs
            for cidx in group:
                _, idx_list = tet_clusters_sorted[cidx]
                idx_list = list(idx_list)
                pos = positions_tet[idx_list]
                long_base = pos[:,0]*row_dir[0] + pos[:,1]*row_dir[1]
                row_long_list.append(long_base)
                row_pos_list.append(pos)
                
            for seg_k in range(len(row_long_list)):
                order = np.argsort(row_long_list[seg_k])
                row_long_list[seg_k] = row_long_list[seg_k][order]
                row_pos_list[seg_k]  = row_pos_list[seg_k][order]
                if seg_k == 1: #apply shift to only the second cluster (if there is one)
                    row_long_list[seg_k]=row_long_list[seg_k]+half_diag  #along proj shift
                    row_pos_list[seg_k]=row_pos_list[seg_k]- [0,Ly,0]           #positional shift
            
            # ---- concatenate segments ----
            row_long = np.concatenate(row_long_list)      # shape (N,)
            row_pos  = np.vstack(row_pos_list)           # shape (N,3)

            N = len(row_long)
            if N < 2:
                vacancy_pairs = []
            else:
                # ---- sort the whole row along long ----
                order = np.argsort(row_long)
                long_sorted = row_long[order]
                pos_sorted  = row_pos[order]

                # ---- compute internal gaps ----
                diffs = np.diff(long_sorted)    # gaps between consecutive atoms
               
    
                # ---- detect base spacing ----
                base_spacing = np.median(diffs)
                # ---- wrap gap across boundary ----
                wrap_gap = long_sorted[0] - long_sorted[-1] + L_diag
                # collect all gaps: diffs + wrap_gap (at end)
                all_gaps = list(diffs) + [wrap_gap]
                # the corresponding index pairs (i,j)
                gap_pairs = [(i, i+1) for i in range(N-1)] + [(N-1, 0)]
                # ---- vacancy detection ----
                vacancy_pairs = [ (i, j)
                                for (gap, (i,j)) in zip(all_gaps, gap_pairs)
                                if gap > 1.5 * base_spacing ]


    
       
         

            for (i, j) in vacancy_pairs:
                p_i = pos_sorted[i]
                p_j = pos_sorted[j]
                if j !=0:
                     # midpoint in real 3D space (unwrapped along the row)
                    vac = 0.5 * (p_i + p_j)
                else:
                    #build periodic image of 0
                    vac= 0.5 * (p_i+ p_j + np.array([Lx,Ly,0]))
              
                # now wrap into the periodic cell
                x = vac[0] % Lx
                y = vac[1] % Ly
                z = vac[2]          # z does not need wrapping

                vacancies.append([x,y,z])
                
   
    return (
        z_oct,
        fe_oct_centers,
        z_tet,
        fe_tet_centers,
        vacancies
    )



def fe_rows_analysis_full(
        fe_top_sel, fe_bot_sel,
        fe_tet_top_sel, fe_tet_bot_sel,
        row_cluster_threshold,
        top_row_normal, bottom_row_normal, box
    ):
    """
    Run the Fe row + vacancy analysis for both surfaces.
    """

    # ---- TOP SURFACE ----
    (z_top,
     fe_top_centers,
     z_tet_top,
     fe_tet_top_centers,
     vac_top) = analyze_fe_rows_for_surface(
                    fe_top_sel,
                    fe_tet_top_sel,
                    row_cluster_threshold,
                    top_row_normal, box
                )

    # ---- BOTTOM SURFACE ----
    (z_bot,
     fe_bot_centers,
     z_tet_bot,
     fe_tet_bot_centers,
     vac_bot) = analyze_fe_rows_for_surface(
                    fe_bot_sel,
                    fe_tet_bot_sel,
                    row_cluster_threshold,
                    bottom_row_normal, box
                )

    # return everything cleanly
    return (
        z_top, z_bot,
        fe_top_centers, fe_bot_centers,
        z_tet_top, z_tet_bot,
        fe_tet_top_centers, fe_tet_bot_centers,
        vac_top, vac_bot
    )




# ---------- core analysis ----------
def analyze_with_mda(u, 
                     atom_names, fe_buffer=1, row_cluster_threshold=2, 
                     top_row_normal = [1/sqrt(2), 1/sqrt(2), 0],
                     bottom_row_normal = [1/sqrt(2), -1/sqrt(2), 0],
                     surf_bounds=(-1.0, 1.0), 
                     layer_bounds=(1.0, 3.5),
                     dipole_theta_boundary=[45,160],
                     h_cutoff=1.2,
                     onrow_fraction=0.125,
                     vac_cutoff=3,
                     wat_on_vac_h_cutoff=1.3,
                  ):


    
    f_oxy_top  = open("oxy_table_top.csv", "w")
    f_oxy_bot  = open("oxy_table_bot.csv", "w")

    f_sum_top  = open("oxy_summary_top.csv", "w")
    f_sum_bot  = open("oxy_summary_bot.csv", "w")

    
    f_oxy_top.write("frame,atom_index,region,h_count,h_indx,row_label,species\n")
    f_oxy_bot.write("frame,atom_index,region,h_count,h_indx,row_label,species\n")
    f_sum_top.write("frame,dissociated_OH,0_diss,1_diss,2_diss,remaining_OH,wat_on_vac,dipole_up,dipole_middle,dipole_down,oct_wat,tet_wat,OH_on_vac,H3O_on_vac,degenerate\n")
    f_sum_bot.write("frame,dissociated_OH,0_diss,1_diss,2_diss,remaining_OH,wat_on_vac,dipole_up,dipole_middle,dipole_down,oct_wat,tet_wat,OH_on_vac,H3O_on_vac,degenerate\n")

   


    for ts in tqdm(u.trajectory, desc="Processing frames"):

        h_atoms = u.select_atoms(f"name {atom_names[0]}")
        o_atoms = u.select_atoms(f"name {atom_names[1]}")
        fe_atoms = u.select_atoms(f"name {atom_names[2]}")
        
        if ts.frame % 10 == 0:
            # --- surface iron selection ---
            (fe_top_sel, fe_bot_sel, z_top, z_bot,
                fe_tet_top_sel, fe_tet_bot_sel, z_tet_top, z_tet_bot) = locate_fe_layers(fe_atoms, fe_buffer)

            (z_top, z_bot,
            fe_top_centers, fe_bot_centers,
            z_tet_top, z_tet_bot,
            fe_tet_top_centers, fe_tet_bot_centers,
            vac_top, vac_bot) = fe_rows_analysis_full(
            fe_top_sel, fe_bot_sel,
            fe_tet_top_sel, fe_tet_bot_sel,
            row_cluster_threshold,
            top_row_normal, bottom_row_normal, box=ts.dimensions
            )
            
   


        classify_surface_oxygens(
            h_atoms,
            o_atoms,
            z_top=z_top,
            z_bot=z_bot,
            surf_bounds=surf_bounds,
            layer_bounds=layer_bounds,
            frame=ts.frame,
            file_top=f_oxy_top,
            file_bot=f_oxy_bot,
            f_sum_top=f_sum_top,
            f_sum_bot=f_sum_bot,
            dipole_theta_boundary=dipole_theta_boundary,
            h_cutoff=h_cutoff,
            fe_top_centers=fe_top_centers,
            fe_bot_centers=fe_bot_centers,
            onrow_fraction=onrow_fraction,
            vac_cutoff=vac_cutoff,
            wat_on_vac_h_cutoff=wat_on_vac_h_cutoff,
            vac_top=vac_top,
            vac_bot=vac_bot,
            box=ts.dimensions
        )



    f_oxy_top.close()  
    f_oxy_bot.close() 

    f_sum_top.close() 
    f_sum_bot.close()



def main():

    parser = argparse.ArgumentParser(description="Label O atoms near slab surfaces: surface-first_layer distinction, H-count, oct/tet analysis.")
    parser.add_argument("dumpfile", help="LAMMPS dump file (coordinates only)")
    parser.add_argument("--fe-buffer", type=float, default=1,
                        help="Fe-surface buffer thickness in Å (default: 1)")
    parser.add_argument("--surface-window", nargs=2, type=float, default=[-1.0, 1.0],
                        metavar=("ZMIN", "ZMAX"),
                        help="dz range (Å) for surface oxygens")
    parser.add_argument("--first-layer-window", nargs=2, type=float, default=[0.5, 3.8],
                        metavar=("ZMIN", "ZMAX"),
                        help="dz range (Å) for first water layer oxygens")
    parser.add_argument("--dipole-theta-boundary", nargs=2, type=float, default=[45, 145],
                        help="boundary to distinguih up-middle-down theta orientation of dipole of water on a vacancy in degree ")
    parser.add_argument("--h-cutoff", type=float, default=1.2,
                        help="Distance cutoff (Å) for O–H bonding")
    parser.add_argument("--water-on-vacancy-h-cutoff", type=float, default=1.2,
                        help="Distance cutoff (Å) for O–H bonding of the water above dissociation (longer OH)")
    parser.add_argument("--rowgap", type=float, default=2.0,
                    help="clustering row gap in Å (default: 2.0)")  
    parser.add_argument("--frac", type=float, default=0.125,
                    help="fraction of inter-row distance used for oct/tet split (default: 0.125)")
    parser.add_argument("--vac_cutoff", type=int, default=5,
                    help="max distance between tet water above dissociation site (vacancy) and the site (default: 5 AA)")
    
    #parser args
    args = parser.parse_args()
    dumpfile = args.dumpfile
    fe_buffer = args.fe_buffer
    surf_bounds=tuple(args.surface_window)
    layer_bounds=tuple(args.first_layer_window)
    h_cutoff=args.h_cutoff
    row_cluster_threshold=args.rowgap
    onrow_fraction=args.frac
    wat_on_vac_h_cutoff=args.water_on_vacancy_h_cutoff
    vac_cutoff=args.vac_cutoff
    dipole_theta_boundary=args.dipole_theta_boundary
    print("surf_bounds:", surf_bounds)
    print("layer_bounds:", layer_bounds)
    print("hydrogen cutoff:", h_cutoff)
    print("water on vacancy hydrogen cutoff:", wat_on_vac_h_cutoff)
    print("water on vacancy dipole boundary:", dipole_theta_boundary)
    print(f"Using Fe buffer = {fe_buffer:.3f} Å")
    print(f"Using Fe cluster threshold = {row_cluster_threshold:.3f} Å")
    print(f"Using vacancy threshold = {vac_cutoff:.3f} Å")
    print(f"Loading dump: {dumpfile}")
    
    atom_names= ["H", "O", "Fe"]
    u = build_topo(args.dumpfile,  atom_names)

    print(f"Universe ready: {len(u.atoms)} atoms, {u.trajectory.n_frames} frames")

    analyze_with_mda(u, 
                     atom_names=atom_names,
                     fe_buffer=fe_buffer, 
                     row_cluster_threshold=row_cluster_threshold,
                     surf_bounds=surf_bounds,
                     dipole_theta_boundary=dipole_theta_boundary,
                     layer_bounds=layer_bounds,
                     h_cutoff=h_cutoff,
                     onrow_fraction=onrow_fraction,
                     vac_cutoff=vac_cutoff,
                     wat_on_vac_h_cutoff=wat_on_vac_h_cutoff
                     )

    print("\nSaved:")
    #print("  surface_row_centers.csv")
    #print("  surface_summary.csv")
    print("  oxy_table_top.csv")
    print("  oxy_table_bot.csv")
    print("  oxy_summary_top.csv")
    print("  oxy_summary_bot.csv")
if __name__ == "__main__":
    main()
