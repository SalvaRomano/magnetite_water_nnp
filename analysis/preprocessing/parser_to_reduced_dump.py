# %%
import numpy as np
import matplotlib.pyplot as plt
import re
import h5py
import ase.io as aio



# %%
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

import os


def  make_folder(folder_name):
     if not os.path.exists(folder_name):
        os.mkdir(folder_name)
        
def make_empty_folder(folder_name):
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
    else:
        try:
            # List all files in the directory
            files = os.listdir(folder_name)

            # Iterate over each file and remove it
            for file in files:
                file_path = os.path.join(folder_name, file)
                if os.path.isfile(file_path):
                    os.remove(file_path)
                    #print(f"File '{file_path}' removed successfully.")
        except OSError as e:
            print(f"Error: {e.strerror}")
               
def create_file(filename):
    with open(filename, 'w'):
        pass 

def find_files_by_name(folder, file_name):
    found_files = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file == file_name:
                found_files.append(os.path.join(root, file))
    
    print(list(np.sort(found_files)))
    return list(np.sort(found_files))



def find_files(start_dir, file_name):
    results = []
    for dirpath, _, filenames in os.walk(start_dir):
        if file_name in filenames:
            results.append(os.path.join(dirpath, file_name))
    print(results)    
    return results

def filter_paths(paths, target_folder):
    return [path for path in paths if target_folder in os.path.normpath(path).split(os.path.sep)]

# %%

def wrap_boundary_xyz(H, L_x, L_y, L_z):
    
    L=[L_x,L_y, L_z]
    string_coord=['x', 'y']
    
    for n in range(len(H)):
        for i in range(len(H[n])):
            if  H[n][i]<0:
                print(f"{string_coord[i]} out of 0, {H[n][i]} wrapped in {L[i]+H[n][i]}" )
                H[n][i] = L[i]+H[n][i]
                
            elif H[n][i]>=L[i]:
                print(f"{string_coord[i]} out of L_{string_coord[i]}, {H[n][i]} wrapped in {-L[i]+H[n][i]}" )
                H[n][i] = -L[i]+H[n][i]
                
def wrap_boundary_xy(H, L_x, L_y):
    
    L=[L_x,L_y]
    string_coord=['x', 'y']
    
    for n in range(len(H)):
        for i in range(len(H[n])-1):
            if  H[n][i]<0:
                print(f"{string_coord[i]} out of 0, {H[n][i]} wrapped in {L[i]+H[n][i]}" )
                H[n][i] = L[i]+H[n][i]
                
            elif H[n][i]>=L[i]:
                print(f"{string_coord[i]} out of L_{string_coord[i]}, {H[n][i]} wrapped in {-L[i]+H[n][i]}" )
                H[n][i] = -L[i]+H[n][i]

# %%
def name_reduced_dump_files(dump_file_folder, reduction_mode):
    if reduction_mode == 'O':
        file_u=dump_file_folder+'./O_u.dump'
        file_d=dump_file_folder+'./O_d.dump'
    elif reduction_mode == 'W':
        file_u=dump_file_folder+'./W_u.dump'
        file_d=dump_file_folder+'./W_d.dump'
    elif reduction_mode == 'W+Fe':
        file_u=dump_file_folder+'./W+Fe_u.dump'
        file_d=dump_file_folder+'./W+Fe_d.dump'
    else:
        print("ERROR: reduction_mode must be 'O', 'W' or 'W+Fe'")
    
    return file_u, file_d



def write_dump_header(file, timestep, N,  L_x, L_y , L_z, wrap_bool):
    #print("writing header")
    file.write("ITEM: TIMESTEP \n")   
    file.write(f"{timestep} \n")   
    file.write("ITEM: NUMBER OF ATOMS \n")  
    file.write(f"{N} \n")
    file.write("ITEM: BOX BOUNDS pp pp pp \n")  
    file.write( f"0.0 {L_x} \n" )
    file.write( f"0.0 {L_y} \n" )
    file.write( f"0.0 {L_z} \n" )
    if wrap_bool:
        file.write("ITEM: ATOMS id type x y z \n")  
    else:  
        file.write("ITEM: ATOMS id type xu yu zu \n")
        
        
#all position = [H, O, Fe]
def write_slice_to_file( file, all_positions, timestep, L_x, L_y , L_z, wrap_bool):
    
    flattend_positions=[]
    
    #N_species=len(all_positions)
    #print("N_types =" , N_species)
    N_atom_per_species=[]
    for all_pos in all_positions:
        #print(all_pos)
        N_atom_per_species.append(len(all_pos))
        for pos in all_pos:
            flattend_positions.append(pos)
    #print("N_atom_per_species =", N_atom_per_species)  
    N_tot=np.sum(N_atom_per_species)
    #print("N_tot =", N_tot)
    
    with open(file, 'a') as append_file:
        #write header
        write_dump_header(append_file, timestep, N_tot,  L_x, L_y , L_z, wrap_bool)
        
        for index, new_line in enumerate(flattend_positions):
            # convert numbers to strings and join them into one line
            string_line = ' '.join(map(str, new_line) ) + '\n' 
            
            if index < N_atom_per_species[0]:
                string_line = f" 1 {string_line}"
            elif index < N_atom_per_species[0]+N_atom_per_species[1]:
                string_line = f" 2 {string_line}"
            else:
                string_line = f" 3 {string_line}"
            string_line = f"{index + 1} {string_line}"
                             
            append_file.write(string_line)
            
            
def clear_variables():
    return [],  [], [], [], [], [], [], [], [], []

# %%
def all_reduced_positions( H_u ,O_u, OS_u , FeS_u ,FeS2_u, reduction_mode):
    if reduction_mode == 'O':
       all_positions=[[],O_u,[]]
    elif reduction_mode == 'W':
        all_positions=[H_u, O_u,[]]
    elif reduction_mode == 'W+Fe':
        all_positions=[H_u, np.concatenate((O_u,OS_u),axis=0) , np.concatenate((FeS_u ,FeS2_u),axis=0)]  
    return  all_positions       
                 

# %%

                                            
def h5py_write_reduced_dump(path_h5py, surface_selection,  dump_file_folder,  reduction_mode , L_x, L_y , L_z):
    wrap_bool=True


    file_u, file_d = name_reduced_dump_files(dump_file_folder, reduction_mode)
     
    if surface_selection =="u":
        file=file_u
    else:
        file=file_d
        
    create_file(file)

    t=0
    h5py_str='.h5py'
    with h5py.File(path_h5py+"H_"+surface_selection+h5py_str, 'r') as f_H, h5py.File(path_h5py+"O_"+surface_selection+h5py_str, 'r') as f_O, h5py.File(path_h5py+"OS_"+surface_selection+h5py_str, 'r') as f_OS, h5py.File(path_h5py+"FeS_"+surface_selection+h5py_str, 'r') as f_FeS, h5py.File(path_h5py+"FeS2_"+surface_selection+h5py_str, 'r') as f_FeS2:
        for k in f_H.keys(): #main time loop

            H = f_H[k][:]
            O = f_O[k][:]
            OS = f_OS[k][:]   
            FeS=f_FeS[k][:]   
            FeS2=f_FeS2[k][:]  
            #WRITE BLOCK
            #select which atom to print
            all_positions=all_reduced_positions( H ,O, OS , FeS ,FeS2, reduction_mode)
           
            #write time frame to file
            write_slice_to_file( file, all_positions, t, L_x, L_y , L_z, wrap_bool)
            t+=1

        

# %%

                                            
def write_reduced_dump(file_trj, equi_time, dump_file_folder, wrap_bool, periodic_image_bool, reduction_mode , L_x, L_y , L_z, N_water,  N_H, N_O, N_Fe,  N_OS, N_FeS, N_FeS2):

    pattern= r"TIMESTEP"
    real_timestep=0
    
    line_count=0
    header=9

    file_u, file_d = name_reduced_dump_files(dump_file_folder, reduction_mode)
     
    create_file(file_u)
    create_file(file_d)

    H_u, H_d ,O_u, OS_u ,OS_d , O_d ,FeS_u ,FeS2_u, FeS2_d ,FeS_d= [],  [], [], [], [], [], [], [], [], []
    # Read the trj file. set the header and the pattern of starting of a config
    t=0
    with open(file_trj, 'r') as file:

        for count, line in enumerate(file):
            #NEW STRUCTURE
           
            
            if  re.search(pattern, line):
                t+=1
                #WRITE BLOCK
                if t >= equi_time:
                    #select which atom to print
                   
                    all_positions_u=all_reduced_positions( H_u ,O_u, OS_u , FeS_u ,FeS2_u, reduction_mode)
                    all_positions_d=all_reduced_positions( H_d ,O_d, OS_d , FeS_d ,FeS2_d, reduction_mode)
                    #write time frame to file
                    write_slice_to_file( file_u, all_positions_u, real_timestep, L_x, L_y , L_z, wrap_bool)
                    write_slice_to_file( file_d, all_positions_d, real_timestep, L_x, L_y , L_z, wrap_bool)
                #CLEAR
                line_count=0
                initial_struct_line=count
                H_u, H_d ,O_u, OS_u ,OS_d , O_d ,FeS_u ,FeS2_u, FeS2_d ,FeS_d=clear_variables()

            #READING TRJ DUMP BLOCK   (LOW COVERAGE)
            #timestep
            if count == initial_struct_line+1:
                real_timestep=line.split()[0]
                #print(real_timestep)
            #position block
            if line_count >= header:           
                split_line=line.split()
                #print(line_count)
                if wrap_bool:
                    position= [float(split_line[2]),float(split_line[3]),float(split_line[4])]
                else:
                    if periodic_image_bool:  #PERIODIC IMAGE 
                        #print('periodic image convenction: unwrap traj is x_unwr=x+Lx*ix')
                 
                        position= [float(split_line[2])+  L_x*int(split_line[5]),   float(split_line[3])  +  L_y*int(split_line[6])  ,float(split_line[4])]
                        #print(position)
                    else:
                        position= [float(split_line[5]),float(split_line[6]),float(split_line[4])]

                #Hydrogen up and down
                if line_count<header+N_H/2:
                    H_u.append(position)
                elif line_count<header+N_H:
                    H_d.append(position)

                #Oxygen water up, surface up, (skip bulk) surface down and water down  
                elif line_count<header+N_H+N_water/2:
                    O_u.append(position)
                elif line_count<header+N_H+N_water/2+N_OS:
                    OS_u.append(position)
                elif line_count<header+N_H+ N_O -N_water/2 -N_OS: 
                    pass 
                elif line_count<header+N_H+ N_O -N_water/2 :
                    OS_d.append(position)
                elif line_count<header+N_H+ N_O:
                    O_d.append(position)

                #Fe surface up, second layer up, (skip bulk); down second layer down
                elif line_count<header+N_H+N_O+N_FeS:
                    FeS_u.append(position)
                elif line_count<header+N_H+N_O+N_FeS+N_FeS2:
                    FeS2_u.append(position)
                elif line_count<header+N_H+N_O+N_Fe - N_FeS-N_FeS2 :
                    pass
                elif line_count<header+N_H+N_O+N_Fe -N_FeS:
                    FeS2_d.append(position)
                elif line_count<header+N_H+N_O+N_Fe:
                    FeS_d.append(position)

            line_count+=1
        

# %%

# RUN IN PLACE

#path_run='/data/home/salvatorer26/projects/ferrite_oxide/production/production_runs/replicated_1_2'


dump_file_folder="./reduced_dump/"
file_trj='./trjmagnetite.dump'





wrap_bool=False         # keep false to unwrap
periodic_image_bool=True  #(if wrap_bool=False ) True if the file format is x y z ix iy (False if x y z xy yu) z always fixed
equi_time=2           #equilibration time (not printing initial N steps)
time_digit_lenght=9    #maximal digit of time
reduction_mode='W'    #O for only oxygens of water, W for water (O and H at surface), W+Fe  for also first layer of irons




print("Working in ", file_trj)
f, L_x, L_y, L_z, diag_lenght, N_water, N_H, N_O, N_Fe, N_OS, N_FeS, N_FeS2= get_initial_info(file_trj)
#ATTENTION MAKE EMPTY THE FOLDER
make_empty_folder(dump_file_folder)
write_reduced_dump(file_trj,equi_time,  dump_file_folder, wrap_bool, periodic_image_bool, reduction_mode , L_x, L_y , L_z, N_water,  N_H, N_O, N_Fe,  N_OS, N_FeS, N_FeS2)


