# %%
import numpy as np
import matplotlib.pyplot as plt
import re
import os
import ase.io as aio

# %%
#SELECT PARAMETERS

path_of_runs='./'

file_trj=path_of_runs + 'trjmagnetite.dump'             #original trj file (lammps)
h5py_files_folder=path_of_runs + "trajectory_h5py_files/" #process trj file
folder_hist=path_of_runs + "histogram_bond_angular/"       #saving folder


# Define binning
n_bins=100  #common number of bins (for OH bond length r and theta)

#OH length binning
rho_o=0.9    #OH bond bins initial
delta=1.3    #OH bond bins final


#select how many water layers are defined after the first two
N_extra_layer=14  #equally spaced after first two layers (7.2 A) (full == 7 ,  full_big == 10, full_huge == 14 )






# %%

def  make_folder(folder_name):
     if not os.path.exists(folder_name):
        os.mkdir(folder_name)
        
        
        
def create_file(filename):
    with open(filename, 'w'):
        pass 


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

def find_folder_by_name(root_path, folder_name):
    folders=[]
    
    for root, dirs, files in os.walk(root_path):
        if folder_name in dirs:
             folders.append(os.path.join(root, folder_name))
    return np.sort(folders)


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

#takes: atoms in a layer
#determine exactly the surfaces positions with average of z of all atom in the surface: 
def get_layer_z(A):
    z=A[:,2]
    return np.mean(z), max(z), min(z)

def slice_height(A, z_max):
    tmp=[]
    for a in A:
        if a[2]<=z_max:
            tmp.append(a)
    return np.array(tmp)


# %%
def vector_distance_pbc(u,v , L_x, L_y, L_z):
    difference=[]
    
    dx=u[0]-v[0]
    dy=u[1]-v[1]
    dz=u[2]-v[2]
    
    if np.abs(dx) <= L_x/2:
        difference.append(dx)
    elif dx>0:
        difference.append(dx - L_x)
    else:
        difference.append(dx + L_x)
    
    if np.abs(dy) <= L_y/2:
        difference.append(dy)
    elif dy>0:
        difference.append(dy- L_y)
    else:
        difference.append(dy + L_y)
    
    if np.abs(dz) <= L_z/2:
        difference.append(dz)
    elif dz>0:
        difference.append(dz- L_z)
    else:
        difference.append(dz+ L_z)
    
    return np.array(difference)


import math
def cartesian_internal(u, phi_pi_periodicity ):  
    x=u[0]
    y=u[1]
    z=u[2]
    # Calculate radial distance
    r = math.sqrt(x**2 + y**2 + z**2)

    # Calculate polar angle (theta)
    theta = math.acos(z / r) if r != 0 else 0

    # Calculate azimuthal angle (phi)
    phi = math.atan2(y, x)
    if phi_pi_periodicity:
        if phi<0:
            phi = phi+np.pi

    return np.array([r, theta, phi])   

# %%
def get_O_cartesian_OH_internal(H_u_indx, O, H, L_x,L_y, z_max_vacuum, t):
    water_dipole=[]
    OHHH=[]
    OH_water=[]
    OH=[]
    O_alone=[]
    
    
    for i , H_index in enumerate(H_u_indx):

        if len(H_index)==0: 
                #if  O[i][2]>3.7:
                #       print(t,"attention! z O alone =",  O[i][2] )
                O_alone.append(O[i])

        elif len(H_index)==1: 

                OH_internal=cartesian_internal(vector_distance_pbc(H[H_index[0]],O[i]  , L_x, L_y, z_max_vacuum))
                OH.append([O[i],OH_internal])
                
                #if  O[i][2]>3.7:
                #       print(t,"attention! z OH =",  O[i][2] )
                
        elif len(H_index)==2:
                OH_cartesian_0=vector_distance_pbc(H[H_index[0]], O[i] , L_x, L_y, z_max_vacuum)
                OH_cartesian_1=vector_distance_pbc(H[H_index[1]], O[i] , L_x, L_y, z_max_vacuum)
                dipole=OH_cartesian_0+OH_cartesian_1
                
                OH_internal_0=cartesian_internal(OH_cartesian_0)
                OH_internal_1=cartesian_internal(OH_cartesian_1)
                dipole_internal=cartesian_internal(dipole)
                
                OH_water.append([O[i],OH_internal_0])
                OH_water.append([O[i],OH_internal_1])
                water_dipole.append([O[i],dipole_internal] )

        elif len(H_index)==3:
                #if  O[i][2]>3.7:
                #       print(t,"attention! z OHHH =",  O[i][2] )
                OH_cartesian_0=vector_distance_pbc(H[H_index[0]], O[i] , L_x, L_y, z_max_vacuum)
                OH_cartesian_1=vector_distance_pbc(H[H_index[1]], O[i] , L_x, L_y, z_max_vacuum)
                OH_cartesian_2=vector_distance_pbc(H[H_index[2]], O[i] , L_x, L_y, z_max_vacuum)
                OH_internal_0=cartesian_internal(OH_cartesian_0)
                OH_internal_1=cartesian_internal(OH_cartesian_1)
                OH_internal_2=cartesian_internal(OH_cartesian_2)
                
                OHHH.append([O[i],OH_internal_0])
                OHHH.append([O[i],OH_internal_1])
                OHHH.append([O[i],OH_internal_2])
                
        else:
              print( "attention! found H*O=", H_index, "at z=", O[i][2])  
    return water_dipole , OH_water , OH, OHHH, O_alone

# %%


def get_O_cartesian_OH_internal(H_u_indx, O, H, L_x,L_y, z_max_vacuum, t, phi_pi_periodicity):
    water_dipole=[]
    OHHH=[]
    OH_water=[]
    OH=[]
    O_alone=[]


    for i , H_index in enumerate(H_u_indx):

        if len(H_index)==0:
                #if  O[i][2]>3.7:
                #       print(t,"attention! z O alone =",  O[i][2] )
                O_alone.append(O[i])

        elif len(H_index)==1:

                OH_internal=cartesian_internal(vector_distance_pbc(H[H_index[0]],O[i]  , L_x, L_y, z_max_vacuum), phi_pi_periodicity)
                OH.append([O[i],OH_internal])

                #if  O[i][2]>3.7:
                #       print(t,"attention! z OH =",  O[i][2] )

        elif len(H_index)==2:
                OH_cartesian_0=vector_distance_pbc(H[H_index[0]], O[i] , L_x, L_y, z_max_vacuum)
                OH_cartesian_1=vector_distance_pbc(H[H_index[1]], O[i] , L_x, L_y, z_max_vacuum)
                dipole=OH_cartesian_0+OH_cartesian_1

                OH_internal_0=cartesian_internal(OH_cartesian_0, phi_pi_periodicity)
                OH_internal_1=cartesian_internal(OH_cartesian_1, phi_pi_periodicity)
                dipole_internal=cartesian_internal(dipole, phi_pi_periodicity)

                OH_water.append([O[i],OH_internal_0])
                OH_water.append([O[i],OH_internal_1])
                water_dipole.append([O[i],dipole_internal] )

        elif len(H_index)==3:
                #if  O[i][2]>3.7:
                #       print(t,"attention! z OHHH =",  O[i][2] )
                OH_cartesian_0=vector_distance_pbc(H[H_index[0]], O[i] , L_x, L_y, z_max_vacuum)
                OH_cartesian_1=vector_distance_pbc(H[H_index[1]], O[i] , L_x, L_y, z_max_vacuum)
                OH_cartesian_2=vector_distance_pbc(H[H_index[2]], O[i] , L_x, L_y, z_max_vacuum)
                OH_internal_0=cartesian_internal(OH_cartesian_0, phi_pi_periodicity)
                OH_internal_1=cartesian_internal(OH_cartesian_1, phi_pi_periodicity)
                OH_internal_2=cartesian_internal(OH_cartesian_2, phi_pi_periodicity)

                OHHH.append([O[i],OH_internal_0])
                OHHH.append([O[i],OH_internal_1])
                OHHH.append([O[i],OH_internal_2])

        else:
              print( "attention! found H(>4)O=", H_index, "at z=", O[i][2])
    return water_dipole , OH_water , OH, OHHH, O_alone






def get_O_cartesian_OH_internal_water_plus_H3O(H_u_indx, O, H, L_x,L_y, z_max_vacuum, t, phi_pi_periodicity):
    water_dipole=[]
 
    OH_water=[]
    OH=[]
    O_alone=[]
    
    
    for i , H_index in enumerate(H_u_indx):

        if len(H_index)==0: 
                #if  O[i][2]>3.7:
                #       print(t,"attention! z O alone =",  O[i][2] )
                O_alone.append(O[i])

        elif len(H_index)==1: 

                OH_internal=cartesian_internal(vector_distance_pbc(H[H_index[0]],O[i]  , L_x, L_y, z_max_vacuum), phi_pi_periodicity)
                OH_water.append([O[i],OH_internal])
                
                #if  O[i][2]>3.7:
                #       print(t,"attention! z OH =",  O[i][2] )
                
        elif len(H_index)==2:
                OH_cartesian_0=vector_distance_pbc(H[H_index[0]], O[i] , L_x, L_y, z_max_vacuum)
                OH_cartesian_1=vector_distance_pbc(H[H_index[1]], O[i] , L_x, L_y, z_max_vacuum)
       
                
                OH_internal_0=cartesian_internal(OH_cartesian_0, phi_pi_periodicity)
                OH_internal_1=cartesian_internal(OH_cartesian_1, phi_pi_periodicity)

                
                OH_water.append([O[i],OH_internal_0])
                OH_water.append([O[i],OH_internal_1])
  
                
        elif len(H_index) >= 3:
                #if  O[i][2]>3.7:
                #       print(t,"attention! z of H", len(H_index), "H =",  O[i][2] )
                
                
                for h in range(len(H_index)):
                
                        OH_cartesian=vector_distance_pbc(H[H_index[h]], O[i] , L_x, L_y, z_max_vacuum)
                        OH_internal=cartesian_internal(OH_cartesian, phi_pi_periodicity)
                        OH_water.append([O[i],OH_internal])
     
    return water_dipole , OH_water , OH,  O_alone

# %%
def OH_z_layer_grouping(O_bond_format,  layer_bins):
    

    # Group arrays based on bins
    grouped_arrays = {category: [] for category in layer_bins}

    for combo in O_bond_format:
        O_z=combo[0][2]
        OH_internal=combo[1]
        for category, condition in layer_bins.items():
            if condition(O_z):
                grouped_arrays[category].append(OH_internal)
                break  # Break out of the loop after categorizing the array

    return grouped_arrays


# %%

def dict_bins_maker(names,layers):  #name of the z layers and layers contain  min and max value
    layer_bins={}
    

    for i, key in enumerate(names):
        layer_bins[key] = (lambda i=i: lambda x: layers[i] < x <= layers[i + 1])()
    
    return  layer_bins

def init_bins_rho_theta_phi(n_bin, delta, rho_o, phi_pi_periodicity):
    
    rho_bins   = np.linspace(rho_o, delta, n_bin) # arbitrary interval (delta is the cutoff)
    theta_bins = np.linspace(-1, 1, n_bin) 
    if phi_pi_periodicity:
        phi_bins   = np.linspace(0 , np.pi , n_bin)
    else:
        phi_bins   = np.linspace(-np.pi, np.pi , n_bin) 
    
    rho_bins_centers=(rho_bins[1:]+rho_bins[:-1])/2
    theta_bins_centers=(theta_bins[1:]+theta_bins[:-1])/2
    phi_bins_centers=(phi_bins[1:]+phi_bins[:-1])/2   
    
    return   rho_bins, theta_bins, phi_bins , rho_bins_centers , theta_bins_centers, phi_bins_centers


def init_combo_hist(layer_bins, n_bin):
    #histograms
   
    combo_hist  = {}

    
    # lenght to normalize
    
    hist_normalization= {}

    
    for layer in layer_bins.keys():
      
        combo_hist[layer]  =[[0]*(n_bin-1)]*(n_bin-1)

        hist_normalization[layer] = 0

    return combo_hist, hist_normalization



# %%

def sum_combo_hist(OH_grouped, combo_hist, theta_bins, rho_bins , hist_normalization,  density_bool ):

    for layer, bonds in OH_grouped.items():
        if len(bonds) !=0:
            theta=np.transpose(np.array(bonds))[1]
            rho=np.transpose(np.array(bonds))[0]    
            hist_normalization[layer]+=len(theta)
            combo_hist[layer]+=np.histogram2d( np.cos(theta),  rho,    bins=[ theta_bins,   rho_bins] , density=density_bool )[0]
          
    return combo_hist, hist_normalization

# %%



def sum_hist_theta_phi(OH_grouped, theta_hist, phi_hist, theta_bins, phi_bins ,density_bool,phi_sin_weight_bool ):

    for layer, bonds in OH_grouped.items():
        if len(bonds) !=0:
            theta=np.transpose(np.array(bonds))[1]
            phi=np.transpose(np.array(bonds))[2]    
            
            theta_hist[layer]+=np.histogram(np.cos(theta), bins=theta_bins,  density=density_bool)[0]
            
            
            
            if phi_sin_weight_bool:
                phi_weights=np.sin(theta)
                phi_hist[layer]  +=np.histogram(phi, bins=phi_bins, weights=phi_weights,  density=density_bool )[0]   
            else:
                phi_hist[layer]  +=np.histogram(phi, bins=phi_bins,  density=density_bool )[0]   
                
    return theta_hist, phi_hist   
    

# %%
def safe_concatenate(a,b):
    if len(a)==0:
        return b
    elif len(b)==0:
        return a
    else:
        return np.concatenate((a,b),axis=0) 

# %%

import h5py
from scipy.spatial import cKDTree

def main_bond_length(h5py_files_folder, folder_hist, n_bins,  surface_selection, L_x, L_y, z_max , rho_o ,  delta,   layers, names,   density_bool, phi_pi_periodicity):
    print('working with ',h5py_files_folder )
    print('surface selected =', surface_selection )
    z_max_H=z_max+delta           #max distance from the surface for H (delta above to include the missing H at the cut)
    z_max_vacuum= z_max+2*delta   #another delta of vacuum distance on top to avoid wrong PBC count


    #H tree box for bonding
    box_size = np.array([L_x, L_y, z_max_vacuum])  
    
    

    # Init bins as dictionary with key as layers (all, oct, tet)
    layer_bins = dict_bins_maker(names,layers)
    #initialize bins of rho, theta, phi  and respective bin centers
    rho_bins, theta_bins, phi_bins , rho_bins_centers , theta_bins_centers, phi_bins_centers=init_bins_rho_theta_phi(n_bins, delta, rho_o, phi_pi_periodicity)

    d_rho=(rho_bins[-1]-rho_bins[0])/(n_bins-1)
    d_theta=(theta_bins[-1]-theta_bins[0])/(n_bins-1)
    d_phi=(phi_bins[-1]-phi_bins[0])/(n_bins-1)

    #initialize histograms with zeros
    combo_hist_OH      ,  norm_hist_OH        =    init_combo_hist(layer_bins, n_bins)
    combo_hist_OH_water,  norm_hist_OH_water  =    init_combo_hist(layer_bins, n_bins)
    combo_hist_dipole  ,  norm_hist_dipole    =    init_combo_hist(layer_bins, n_bins)






    #TIME LOOP BEGIN
    t=0
    h5py_str='.h5py'
    with h5py.File(h5py_files_folder+"H_"+surface_selection+h5py_str, 'r') as f_H, h5py.File(h5py_files_folder+"O_"+surface_selection+h5py_str, 'r') as f_O, h5py.File(h5py_files_folder+"OS_"+surface_selection+h5py_str, 'r') as f_OS:
        for k in f_H.keys(): #main time loop

            H = f_H[k][:]
            O = f_O[k][:]
            OS = f_OS[k][:]

            #add oxygen at surface to get the OH bond of dissociation
            all_O=np.concatenate((O,OS))

            # translate to 0 the surfaces (and invert the d)
            O=slice_height(all_O, z_max)
            H=slice_height(H, z_max_H)

            #H index is a list of list of H indices: the 0 index is the Oxygen, the second runs from 0 to the N_H that are within a cutoff delta to that O
            tree = cKDTree(data=H, leafsize=H.shape[0], boxsize= box_size)  #makes a tree of Hs
            H_indx = tree.query_ball_point(O, delta)                        #assign every H to an O within the cutoff and make the list of list H


            #extract from the H_indx the Oxygen in cartesian and the OH(s)in polar:
            #the output is then a list of pair [O_vect, OH]

            #dipole , OH_water , OH, O_HHH, O_alone=get_O_cartesian_OH_internal(H_indx, O, H, L_x,L_y, z_max_vacuum, t)   #not periodic distance for the z direction
            dipole , OH_water , OH, OHHH, O_alone=get_O_cartesian_OH_internal(H_indx, O, H, L_x,L_y, z_max_vacuum, t, phi_pi_periodicity)


            #use the info of the z Oxygen to group the OH by their heights: output is dictonary with layer_bins keys and vlues OH in internal
            dipole_z  = OH_z_layer_grouping(dipole,   layer_bins)
            OH_water_z= OH_z_layer_grouping(OH_water, layer_bins)
            OH_z      = OH_z_layer_grouping(OH,       layer_bins)


            #Sum the histograms
            combo_hist_OH      ,  norm_hist_OH,         = sum_combo_hist(OH_z,       combo_hist_OH,       theta_bins, rho_bins,   norm_hist_OH,       density_bool)
            combo_hist_OH_water,  norm_hist_OH_water,   = sum_combo_hist(OH_water_z, combo_hist_OH_water, theta_bins, rho_bins,   norm_hist_OH_water, density_bool )
            combo_hist_dipole  ,  norm_hist_dipole,     = sum_combo_hist(dipole_z,   combo_hist_dipole,   theta_bins, rho_bins,   norm_hist_dipole,   density_bool)

            t+=1
    


    #save histograms
    print('save combo_hist_*** for surface', surface_selection , 'in', folder_hist)
    
    
    make_folder(folder_hist)
    
                #name                                               #bin center                                 #hist               #normalization   
    np.save(folder_hist+"combo_hist_OH_water_"+ surface_selection, (  rho_bins_centers, theta_bins_centers , combo_hist_OH_water   ,  norm_hist_OH_water ))
    np.save(folder_hist+"combo_hist_OH_"+ surface_selection, (  rho_bins_centers, theta_bins_centers , combo_hist_OH   ,  norm_hist_OH))
    np.save(folder_hist+"combo_hist_dipole_"+ surface_selection, (  rho_bins_centers, theta_bins_centers , combo_hist_dipole  ,  norm_hist_OH_water ))
    
    
    
    
    np.save(folder_hist+"layering", layers)
    np.save(folder_hist+'tot_time.txt', t)



# %%
#RUNNING BLOCK
z_max=np.load(h5py_files_folder+'z_max.npy') -delta

print("z_max=" , z_max)


# Histograms Init 


#COMMON LAYERS
layers =np.array([-0.5, 1, 3.6, 7.2])

#EXTRA
layers=np.concatenate    ( (layers, np.linspace(7.2, z_max ,  N_extra_layer +1)[1:] ))
print(layers)
# Define name of bins 
names    =["l_"+str(i) for i in np.arange(0,len(layers)-1)]


f, L_x, L_y, L_z, diag_lenght, N_water, N_H, N_O, N_Fe, N_OS, N_FeS, N_FeS2= get_initial_info(file_trj)
print(h5py_files_folder)

phi_pi_periodicity='True'     # True: sin(theta) weighted average for phi  (2D plot always False!)
density_bool=False           # density histogram directly 




#MAIN
for surface_selection in ['u', 'd']:
    main_bond_length(h5py_files_folder, folder_hist, n_bins,  surface_selection, L_x, L_y, z_max , rho_o ,  delta,   layers, names,   density_bool, phi_pi_periodicity)




