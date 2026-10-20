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

def get_initial_info(trj):
    pattern= r"TIMESTEP"
    header=8
    N_H=0
    N_O=0
    N_Fe=0
    #SCV structure (1uc)
    N_OS_uc=8    # number oxygen in first outer layer/uc
    N_FeS_uc=4   # number of iron in first outer layer/uc
    N_FeS2_uc=3  # number of iron  second/outer layer/uc
    with open(trj, 'r') as file:

        for count, line in enumerate(file):     
            if  re.search(pattern, line) and count !=0 :
                break
                
            split_line=line.split()   
            if count == 3:
                N_atom=float(split_line[0])
            if count == 5:
                L_x=float(split_line[1])
            if count == 7:
                L_z=float(split_line[1])
            #position block
            if count > header:        
                if split_line[1] == '1':
                   N_H+=1
                if split_line[1] == '2':
                   N_O+=1   
                if split_line[1] == '3':
                   N_Fe+=1   
        L_y=L_x
        
        #replica
        for n in range(1,10):
            if 8 < L_x/n < 9:
                n_replica=n
                
        if N_atom != N_H+N_O+N_Fe:
            print("ATOM LOST")
        
        f=n_replica**2
        N_water=int(N_H/2)
        N_Fe_uc=int(N_Fe/f)
        N_O_uc=int((N_O-N_water)/f)
        w_uc=float(N_water/f/2)
        N_OS=N_OS_uc*f
        N_FeS=N_FeS_uc*f
        N_FeS2= N_FeS2_uc*f
        diag_lenght=np.sqrt(L_x**2+L_y**2)
        print("L_x = L_y = ", L_x)
        print("L_z = ", L_z)
        print("Diagonal= ", diag_lenght)
        print("N_O_uc = ",  N_O_uc)
        print("N_Fe_uc =",  N_Fe_uc)
        print("w_uc =",  w_uc)
        print(f"Size factor f = {f}")
        print(f"N_water = {N_water}")
        print(f"N_H = {N_H}")
        print(f"N_O = {N_O}")
        print(f"N_Fe = {N_Fe}")
        print(f"N_atom = {N_atom}")
        print(f"N_O surface = {N_OS}")
        print(f"N_Fe surface = {N_FeS}")
        print(f"N_Fe surface2 = {N_FeS2}")
        
        return(f, L_x, L_y, L_z, diag_lenght, N_water,  N_H, N_O, N_Fe, N_OS, N_FeS, N_FeS2)
    
# %%

def msd_main(path_dump, surface_selection, atom_types, atom_names ):
    u = MDAnalysis.Universe(path_dump + "W_"+surface_selection+".dump", format="LAMMPSDUMP")


    #TOTAL DIFFUSION: Basis vector MSD analysis

    #results == here msd for each atom will be stored      #type 2 is the oxygen
    results_O = np.zeros((len(u.trajectory), len(u.select_atoms(atom_type))))


    #all the trajectory in an array [t,n,[x,y,z]]
    coord_path_O = np.zeros((len(u.trajectory), len(u.select_atoms(atom_type)), 2))
    coord_path_O = np.zeros((len(u.trajectory), len(u.select_atoms(atom_type)), 2))

    for ts in u.trajectory:                                    #sliced to take only xy 
        coord_path_O[ts.frame] = u.select_atoms(atom_type).positions[:,:2]/L_x  ### #expect normalized coordinates



    #msd calculation
    for n in range(len(u.select_atoms(atom_type))):
        results_O[:, n] = tidynamics.msd(coord_path_O[:,n,:])

    average_msd = results_O.mean(axis=1)
    np.savez(msd_folder + atom_names[count]+'_MSD_xy_'+ surface_selection+'.npz',  msd=average_msd)
    np.savez(msd_folder + atom_names[count]+'_MSD_xy_all_atoms_'+ surface_selection+'.npz',  msd=results_O)
    #DIRECTIONAL DIFFUSION 2D XY
    special_directions = [
        [1/np.sqrt(2),  1/np.sqrt(2)],  # x + y
        [1/np.sqrt(2), -1/np.sqrt(2)]   # x - y
    ]

    if surface_selection == 'u':
        surface_n_strings = [ 'para', 'ortho']
    else:
         surface_n_strings = [ 'para', 'ortho']

    for k in range(len(special_directions)):

        along_direction = special_directions[k] 

        results_O = np.zeros((len(u.trajectory), len(u.select_atoms(atom_type))))
        coord_path_O = np.zeros((len(u.trajectory), len(u.select_atoms(atom_type)), 2))

        for ts in u.trajectory:
            for n in range(len(u.select_atoms(atom_type))):
                coord_path_O[ts.frame, n] = np.dot(u.select_atoms(atom_type).positions[n,:2], along_direction )  / L_x * np.array(along_direction) #expect normalized coordinates


        for n in range(len(u.select_atoms(atom_type))):
            results_O[:, n] = tidynamics.msd(coord_path_O[:,n,:])


        average = results_O.mean(axis=1)

        np.savez(msd_folder + atom_names[count] + '_MSD_xy_'           + surface_n_strings[k] +"_" + surface_selection +  '.npz',  msd=average)
        np.savez(msd_folder + atom_names[count] + '_MSD_xy_all_atoms_' + surface_n_strings[k] +"_" + surface_selection +  '.npz',  msd=results_O)

# %%

####MAIN###############
import re


path="./"


path_dump=path+"reduced_dump/"
#path_h5py=path+"trajectory_h5py_files/"
msd_folder=path+"MSD.npz_files/"

        
make_folder(msd_folder)


file_trj=path+'trjmagnetite.dump'


f, L_x, L_y, L_z, diag_lenght, N_water, N_H, N_O, N_Fe, N_OS, N_FeS, N_FeS2= get_initial_info(file_trj)
# %%

atom_types=["type 1", "type 2"]
atom_names=['H', 'O']


for surface_selection in['u','d']:
    print('selected surface:', surface_selection)
    for count, atom_type in enumerate(atom_types):
        print('selected atom', atom_type, '==',atom_names[count] )
        msd_main(path_dump, surface_selection, atom_types, atom_names )


