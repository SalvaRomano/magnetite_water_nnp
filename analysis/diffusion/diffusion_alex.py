# %%
import MDAnalysis
import tidynamics 
import numpy as np
import MDAnalysis.lib.mdamath as mdamath
import h5py
import os
# %%

def  make_folder(folder_name):
     if not os.path.exists(folder_name):
        os.mkdir(folder_name)
        
# %%
path = "/data/home/salvatorer26/projects/ferrite_oxide/production/production_runs/water_1_uc/run_4/trajectory_dump_files/"
path_h5py="/data/home/salvatorer26/projects/ferrite_oxide/production/production_runs/water_1_uc/run_4/trajectory_h5py_files"
msd_folder=path+"msd_npz_files"

L_x=25.08508
surface_selection='u'

        
make_folder(msd_folder)


# %%






u = MDAnalysis.Universe(path + "W_" +surface_selection +".dump", format="LAMMPSDUMP")

#Slicing position take only x and y
#t =0
#for ts in u.trajectory:
    
    #print(u.select_atoms("type 2").positions[:,:2]/L_x)
    #print(u.select_atoms("type 1").positions[:,:2])
#    print(u.atoms.positions[:,:2]/L_x)
#    t+=1
#    if t > 2:
#        break 


# %%


#TOTAL DIFFUSION: Basis vector MSD analysis




#wresults == here msd for each atom will be stored      #type 2 is the oxygen
results_O = np.zeros((len(u.trajectory), len(u.select_atoms("type 2"))))


#all the trajectory in an array [t,n,[x,y,z]]

coord_path_O = np.zeros((len(u.trajectory), len(u.select_atoms("type 2")), 2))
coord_path_O = np.zeros((len(u.trajectory), len(u.select_atoms("type 2")), 2))
for ts in u.trajectory:
    coord_path_O[ts.frame] = u.select_atoms("type 2").positions[:,:2]/L_x  ### take only x and y 



#msd calculation
for n in range(len(u.select_atoms("type 2"))):
    results_O[:, n] = tidynamics.msd(coord_path_O[:,n,:])

average_msd = results_O.mean(axis=1)
np.savez(path + 'MSD_xyz'+ '.npz',  msd=average_msd)



# %%
#DIRECTIONAL DIFFUSION 

triclinic_vectors =  mdamath.triclinic_vectors(u.dimensions)
#by_hand triclinic_vectors = [[L_x,0,0], ...]

special_directions = [
    triclinic_vectors[0] + triclinic_vectors[1],  # x + y
    triclinic_vectors[0] - triclinic_vectors[1],  # x - y
]

print( triclinic_vectors[0])


# %%
surface_n_strings = [ 'ortho', 'para']

for k in range(len(special_directions)):

    along_direction = special_directions[k] / np.linalg.norm(special_directions[k])  #versor
   

    results_O = np.zeros((len(u.trajectory), len(u.select_atoms("type 2"))))
    coord_path_O = np.zeros((len(u.trajectory), len(u.select_atoms("type 2")), 3))

    for ts in u.trajectory:
        for n in range(len(u.select_atoms("type 2"))):
            coord_path_O[ts.frame, n] = np.dot(u.select_atoms("type 2").positions[n], along_direction) * along_direction / 25.08508 ## Here is ok since z is doesnt matter


    for n in range(len(u.select_atoms("type 2"))):
        results_O[:, n] = tidynamics.msd(coord_path_O[:,n,:])

        
    average = results_O.mean(axis=1)

    np.savez(path + 'MSD_' + surface_n_strings[k] + '.npz',  msd=average)
    np.savez(path + 'MSD_all_atoms' + surface_n_strings[k] + '.npz',  msd=results_O)


### 2D XY

results = np.zeros((len(u.trajectory), len(u.select_atoms("type 2"))))
coord_path_O = np.zeros((len(u.trajectory), len(u.select_atoms("type 2")), 3))

for ts in u.trajectory:
    coord_path_O[ts.frame] = u.select_atoms("type 2").positions / 25.08508  
    coord_path_O[ts.frame, :, 2] = 0

for n in range(len(u.select_atoms("type 2"))):
    results[:, n] = tidynamics.msd(coord_path_O[:,n,:])

average = results.mean(axis=1)
np.savez(path + 'MSD_xy_plane'+ '.npz',  msd=average)

# %%
yo = np.load(path + 'MSD_' + surface_n_strings[0] + '.npz')
yo = np.load(path + 'MSD_all_atomspara.npz')

print(yo["msd"])

import matplotlib.pyplot as plt

#plt.plot(yo["msd"])
yo = np.load(path + 'MSD_ortho.npz')
plt.plot(yo["msd"], color = "red")

yo = np.load(path + 'MSD_para.npz')
plt.plot(yo["msd"], color = "black")

plt.yscale('log')
plt.xscale('log')

# %%
