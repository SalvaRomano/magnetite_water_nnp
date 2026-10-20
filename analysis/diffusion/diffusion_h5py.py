# %%
#import MDAnalysis
import tidynamics 
import numpy as np
#import MDAnalysis.lib.mdamath as mdamath
import h5py

# %%
path = "/data/home/salvatorer26/projects/ferrite_oxide/production/production_runs/replicated_1_2/random_6wuc_1rep/T_300/run_1/trajectory_dump_files/"
path_h5py="/data/home/salvatorer26/projects/ferrite_oxide/production/production_runs/replicated_1_2/random_6wuc_1rep/T_300/run_1/trajectory_h5py_files/"
surface_selection='u'

h5py_str='.h5py'

coord_path_O=[]
with h5py.File(path_h5py+"H_"+surface_selection+h5py_str, 'r') as f_H, h5py.File(path_h5py+"O_"+surface_selection+h5py_str, 'r') as f_O:
    for k in f_H.keys(): #main time loop
      
        
        coord_path_O.append(f_O[k][:])
        
        

# %%














triclinic_vectors =  mdamath.triclinic_vectors(u.dimensions)
#by_hand triclinic_vectors = [[1,0,0], ...]

## self.results = np.zeros((self.n_frames, self.n_particles))

results = np.zeros((len(u.trajectory), len(u.select_atoms("type 2"))))


#self._position_array = np.zeros(
#    (self.n_frames, self.n_particles, self.dim_fac))


coord_path_O = np.zeros((len(u.trajectory), len(u.select_atoms("type 2")), 3))

for ts in u.trajectory:
    coord_path_O[ts.frame] = u.select_atoms("type 2").positions

#for n in range(self.n_particles):
#    self.results.msds_by_particle[:, n] = tidynamics.msd(
#        positions[:, n, :])

print(coord_path_O)

for n in range(len(u.select_atoms("type 2"))):
    results[:, n] = tidynamics.msd(coord_path_O[:,n,:])

average = results.mean(axis=1)

# surface_normals=[
#     np.cross(triclinic_vectors[1], triclinic_vectors[2]),
#     np.cross(triclinic_vectors[0], triclinic_vectors[2]),
#     np.cross(triclinic_vectors[0], triclinic_vectors[1]),
#     #-triclinic_vectors[0] + 2*triclinic_vectors[2],
#     #triclinic_vectors[0]  + 2*triclinic_vectors[2],
# ]

special_directions = [
    triclinic_vectors[0] + triclinic_vectors[1],  # x + y
    triclinic_vectors[0] - triclinic_vectors[1],  # x - y
]

surface_n_strings = ['para', 'ortho']#, 'm102', '102']

for k in range(len(special_directions)):

    along_direction = special_directions[k]
    along_direction = along_direction / np.linalg.norm(along_direction)

    results = np.zeros((len(u.trajectory), len(u.select_atoms("type 2"))))
    coord_path_O = np.zeros((len(u.trajectory), len(u.select_atoms("type 2")), 3))

    for ts in u.trajectory:
        for n in range(len(u.select_atoms("type 2"))):
            coord_path_O[ts.frame, n] = np.dot(u.select_atoms("type 2").positions[n], along_direction) * along_direction


    for n in range(len(u.select_atoms("type 2"))):
        results[:, n] = tidynamics.msd(coord_path_O[:,n,:])

    average = results.mean(axis=1)

    np.savez(path + 'MSD_' + surface_n_strings[k] + '.npz',  msd=average)

# %%
