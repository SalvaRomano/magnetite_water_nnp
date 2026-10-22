MAGNETITE / WATER INTERFACE WITH A COMMITTEE NEURAL-NETWORK POTENTIAL
=====================================================================
Salvatore Romano (University of Vienna)

Fe3O4(001) slabs in contact with water (full and low coverage, pH 11),
simulated with LAMMPS + an n2p2 committee high-dimensional NNP (C-HDNNP)
trained on optB88-vdW + U DFT data with an active-learning loop.

This repository holds inputs, scripts, the trained potential, small processed
results and one snapshot per system. Trajectories (~600 GB) are NOT included;
see "Heavy data" at the bottom.


LAYOUT
------
committee_nnp/          the trained committee (n2p2): input.nn, scaling.data,
                        weights.{001,008,026}.data = member 0, nnp-data-1/ = member 1,
                        input.every10.data + validation.every10.data (every 10th
                        structure of the training / validation sets, for GitHub size limits)
  training/               provenance of the two members (learning curves, logs, test parity)
active_learning/        the committee active-learning loop (see WORKFLOW)
  n2p2/  lammps/  vasp/ (+ mlff/)  input_data/ (+ cleaning/)
  vasp/mlff/              first data generation with VASP on-the-fly ML-MD
md/                     LAMMPS production runs (inputs + last-frame snapshots)
  tools/                  run-tree set-up, restart and simulated-time helpers
structures/             starting configurations (POSCAR / LAMMPS data) of the MD runs
analysis/               trajectory analysis scripts
results/                paper figures


POTENTIAL (committee_nnp)
-------------------------
- Elements H, O, Fe; 2 hidden layers x 25 nodes, tanh; 156 symmetry functions
  (radial + angular), cutoff 12 bohr = 6.35 A, cutoff_type 3.
- Trained 50 epochs, Kalman filter, energies + forces (force_weight 10).
  input.nn contains "normalize_data_set force", but the n2p2 version used ignored it
  ("Unknown keyword ... Data set normalization is not used", nnp-train.log), so no
  normalization constants / output.nn are needed: input.nn + weights are complete.
- The two members are runs 1 and 2 of a 9-run training (VSC, final_committee_training):
    member 0 (weights.*.data)        run_1, epoch 36
    member 1 (nnp-data-1/)           run_2, epoch 49   (md5-identical to the epoch files)
                     E RMSE train/test (meV/atom)   F RMSE train/test (meV/A)
    member 0               0.94 / 0.91                    111 / 104
    member 1               0.94 / 0.86                    111 / 104
  committee_nnp/training/ has learning-curve.out, nnp-train.log, the training input.nn
  (differs only in seed and committee keywords) and test-set parity data for both.
- Training set 26,789 structures, validation 4,807. Reference DFT (VASP): optB88-vdW,
  DFT+U (U=4.5, J=0.89 eV on Fe d), ENCUT 500, spin-polarised ferrimagnetic
  (+-4 muB on Fe), k-points 2x2x1; POTCARs PAW_PBE Fe_pv, O_GW, H_GW.
  Only every 10th structure is in git (2,679 / 481); the full sets are in
  production/production_committee/ on the workstation (see HEAVY DATA).
- Committee: committee_mode validation, "committee_data nnp-data 2".
- Used in LAMMPS with the modified n2p2 build (n2p2-committee-modified) via
    pair_style nnp dir <committee> ... calcd 1000 wricde ... emap "1:H,2:O,3:Fe"
    pair_coeff * * 6.35013
  It writes disagreement.txt when the committee disagreement exceeds the threshold.


WORKFLOW: ONE ACTIVE-LEARNING ITERATION (from the old vademecum)
----------------------------------------------------------------
Each iteration has 4 steps, one folder each in active_learning/.
A modified n2p2 that writes disagreement.txt is required.

1) n2p2 - train the committee and pick the best weights          [n2p2/]   
   - mkdir n2p2/committee_i; copy input.data, validation.data, input.nn there
   - starter.sh            makes run_1..N with different random seeds
   - sbatch job_n2p2_single.sh in each run   (job_series.sh = serial version)
   - learn_weights.ipynb   learning curves -> best epoch
     (get_weights.py / plot_learning_curves.py do the same from the shell)
   - dataset.sh            nnp-dataset with the chosen weights
   - dataset.plus.ipynb    parity plots, writes bad_structures.txt (outliers)
   - input_data/print_config_input_data.ipynb to drop bad structures, retrain if needed
   - nnp-data_maker_for_LAMMPS.sh builds committee_i/ with nnp-data-k/ folders

2) LAMMPS - committee MD and disagreement search                  [lammps/]  
   - md.lmp is the committee-MD template (md.lmp_equi = production template)
   - one folder per starting configuration (data.equi); set nnpDir in md.lmp!
   - sbatch job_committee.sh in each folder
   - get_comm_disag_info.sh   summary of runs and disagreements
   - pruner_disagremeent.sh   (optional) keep disagreements >= 1000 steps apart
   - convert_disagreements_to_poscar.py (or convert_pruned_disag_to_poscar.py)
       -> active_learn_POSCARs/Force_Dis_config_*/POSCAR   (adapt for bulk water!)
   - selector_Force_Dis_config.sh keeps 1 in 10 of them
   - pressure check: log.lammps_to_table.sh + pressure_analysis_LAMMPS.ipynb
       + disagreement_reduce_by_pressure.sh
   - cutter_last.equi.sh      last frame of equi.dump for OVITO

3) VASP - label the new structures                                [vasp/]  
   - copy active_learn_POSCARs/ and setup/ (INCAR, KPOINTS, POTCAR) to VASP/committee_i
   - sbatch job_check.sh (as many times as you want, it picks unfinished folders)
   - get_VASP_committee_info.sh   progress, wrong magnetisation, CGA 250, DAV errors
   - rm_mag_DAV_CGA.sh / touch_bad_mag.sh   remove or flag bad runs
   - OUTCAR_input.data_converter.sh, then cat_input.data_bad_mag.sh -> new input.data
   - select_VASP_test.data_mag_vers.sh (test set);
     pressures: grep -a "external pressure" */OUTCAR

4) input.data - check before training                             [input_data/]  
   - analyis_input.data.ipynb     energies, force norms, lattices; remove bad forces
   - print_config_input_data.ipynb / extract_input_data_from_energy.ipynb
                                  extract structures by index / by energy
   - input.data_to_pdb*.py        movie of a data set for OVITO
   - cleaning/                    duplicates, Fe-force outliers, atom separations
   Then go back to step 1 with the enlarged data set.

test_set_fake_disag_maker.sh (lammps/) writes random "disagreements" to build an
independent test set.


MD SYSTEMS (md)
---------------
LAMMPS 29Oct2020 + n2p2-committee-modified, dt = 0.5 fs, NVE + Langevin
(damp 0.05 ps) on a central slab region only, dump every 500-1000 steps.
Runs are chained through restart_a/b. Each run folder holds md.lmp and the job script
(recovered from the VSC accounts; log.lammps files, 2.3 GB, are not in git);
snapshots/ has the last complete frame of the latest run
(LAMMPS custom dump: id type x y z ...; type 1=H, 2=O, 3=Fe).
Surface cell: (sqrt2 x sqrt2)R45, a = 8.396 A. Slabs are two-sided (up/down).

  folder                        old name                  atoms (H/O/Fe)           time
  full_coverage/small           full_coverage             4352/3200/704, 4x4      ~29 ns (run_1 independent)
  full_coverage/medium          full_coverage_big         6304/4304/800, 4x4      ~29 ns (runs 1-4 independent)
  full_coverage/large           full_coverage_huge        12548/7426/800, 4x4     ~14 ns + 2.5 ns at 350 K
  ph11/small                    ph_11                     full_coverage - 32 H     ~4 ns (equi = NPT)
  ph11/large                    ph_11_huge                20742 atoms              ~3.6 ns
  bulk_truncated/small          bulk_truncated            6304/3792/448            ~3 ns
  bulk_truncated/large          bulk_truncated_big        12514/6897/448           ~1 ns
  slowmo                        slowmo                    = full_coverage/small    1.6 ns, dump every 50 fs
  full_coverage/small/spectra_run  (VSC only)             = full_coverage/small    4 runs, x + v dumped every 4 fs (spectra)
  low_coverage/1x1              replicated_1              1 surface cell           200-500 K, up to ~50 ns
  low_coverage/2x2              replicated_2              2x2 cells                200-400 K
  low_coverage/3x3_and_6x6      replicated_3_6            3x3 (water_1-4), 6x6 (water_5,6,8)   300 K

water_N_uc = N water molecules per surface unit cell, per side.
full_coverage/simulated_times.txt and how_to_concatenate_restarts.txt explain how
restarts were joined (remove overlapping frames, then cat); simulated_times*.txt in
low_coverage/ give the total simulated steps per system (tools/get_production_info.sh).
bulk_truncated/large/ionized: runs of the large bulk-truncated system with ions.
low_coverage/2x2/random_*wuc_2rep: random initial water arrangements (inputs only).
low_coverage/templates/: md.lmp_T_* and job templates used to start the T_* / run_* trees
(tools/starter_T_run.sh); restarts were continued with tools/restarter.sh.
CAVEAT: in low_coverage/{1x1/water_8,9_uc_1rep, 2x2/water_8,9_uc}/T_300 the md.lmp
has temperature_f = 200, so T drifts from 300 K to ~280 K during those runs. The
templates contain both md.lmp_T_300 (300->300) and md.lmp_T_300_200 (300->200 ramp).


ANALYSIS (analysis)
-------------------
Most scripts work inside a run folder that holds trjmagnetite.dump.
Several hard-code the slab geometry (z windows, 64 O/44 Fe vs 72 O/50 Fe per cell)
and old absolute paths: check before running.

preprocessing       parser_to_h5py_single.py: dump -> trajectory_h5py_files/{H,O,OS,FeS,FeS2}_{u,d}.h5py
                    (top/bottom surface split; WARNING: empties the output folder first)
                    parser_to_h5py_ph11.py: same for pH 11 / bulk water
                    parser_to_reduced_dump.py: surface water only -> reduced_dump/W_{u,d}.dump
                    (all of these run inside a run folder; loop over run_* / T_*/run_*)
                    merging_trajectories.ipynb: join reduced dumps of several runs
species_lifetimes   indicator_water_magnetite.py <dump>: per-frame classification of
                    surface O (dissociated/remaining OH, water on vacancy L/R/up/down,
                    H3O, oct/tet water, dipole orientation)
                      -> oxy_table_{top,bot}.csv, oxy_summary_{top,bot}.csv
                      (for i in run_*; do (cd $i && python3 ../indicator_water_magnetite.py trjmagnetite.dump); done)
                    get_all_autocorrelations.py oxy_table_X.csv -> indicator ACFs
                      (run it in run_*/top and run_*/bot); build_indicators.py = helper
                    plot_all_correlation*.ipynb: ACF fits and lifetimes
                    avg_oxy.awk / global_summary.awk: average species counts
structure           density_h5py.ipynb (density profiles, oct/tet split),
                    angular_h5py_file / angular_2D / bond_angular_2D (OH and dipole
                    orientation), plot_hist_h5py_file.ipynb (plots), bond_angular.py,
                    density_h5py_analysis.py (script versions)
diffusion           msd_diffusion_single.py: lateral MSD of O/H, total and along / across
                    the Fe rows -> MSD.npz_files/ (run inside each run folder); plot_save_msd.py
                    (diffusion_h5py.py is broken: imports commented out)
clusters            low_cov_cluster_analysis.ipynb: cluster size and trimer statistics
                    (O-O cutoff 3.36 A) -> cluster_result/
spectra             vibrational SFG spectra from MD (dump -> xyz/xtc, velocity
                    correlation, non-Condon correction, final plots) and VASP
                    polarizability inputs. Needs the external VSFG-magnetite code (not included).
                    dump -> xyz with VSFG-magnetite py/lammps2xyz.py and lammps2xyz_vel.py;
                    average: python3 average_ssVVCF_tree.py --root . --pattern "trj_*/chunk_*/.../ssVVCF_wcf.dat"
density_reweighting    density scans / reruns of bulk water with the NNP
legacy              older ASE-based notebooks


RESULTS (results)
-----------------
Paper figures (png/pdf; inkscape/ = editable svg sources).
magnetite_acs/    low coverage: stable water structures, trimers, densities,
                  dissociation, MSD/diffusion (msd/, diffusion/)
magnetite_jcis/   full coverage: dipole/anisotropy vs film thickness, OH / H-bond
                  analysis, species, pH 11 (ph_11/)
misc/             other figures (anisotropy, dissociation, density, pH 11, theta-rho)
Movies (mp4/gif, ~340 MB) and inkscape files > 20 MB are kept only in ~/Desktop/plots.
The data behind the figures stays in production/production_runs (oxy_*.csv,
histogram_*, MSD.npz_files, data_for_plot, plot_data_anisotropy, cluster_result).


HEAVY DATA (not in git)
-----------------------
On the Vienna workstation (/data/home/salvatorer26/projects/ferrite_oxide/):
  production/production_committee/input.data, validation.data     full training data
  production/production_runs/<old name>/run_*/trjmagnetite.dump   trajectories
  production/production_runs/**/trajectory_h5py_files, histogram_*, MSD.npz_files
  vsc_import/                       files recovered from VSC (training run_1/2 with all
                                    epochs, every log.lammps, VASP_ion_minimization_2)
  production/VASP_ion_minimization_final/                         VASP relaxations (slab, water_6/8;
                                                                  not converged: moments collapse)
  PNAS_script/                                                     spectra data (~2 TB)
Original potential dirs used on VSC: /gpfs/data/fs71576/salvatore/... and
/gpfs/data/fs71796/salvrom/... (same committee as committee_nnp/).
Still only on VSC: all 9 training runs (fs71576/salvatore/ferrite_oxide/final_committee_training),
the spectra_run trajectories (~1.9 TB, fs71796/salvrom/projects/ferrite_oxide/production/
big_systems/production_runs/spectra_run) and the other production trajectories.
