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
                
    return np.sort(found_files)

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
def parser_full_coverage(t_slice,  z_interval_u, z_interval_d,  N_water, N_OS, N_FeS, N_FeS2):
    
    lost_atom_bool=False
    positions=t_slice.get_positions(wrap=True)
    z=positions[:,2]
    #indices species 
    H_indices = np.where(t_slice.get_atomic_numbers()==1)[0]
    O_indices = np.where(t_slice.get_atomic_numbers()==2 )[0]
    Fe_indices = np.where(t_slice.get_atomic_numbers()==3 )[0]
    #indices regions 
    water_indices=np.where( (z > z_interval_u[0] ) |  (z < z_interval_d[0]) )[0]
    
    
    S_u_indices_broad=np.where(  (z < z_interval_u[0] ) &  (z > z_interval_u[1]-0.5) )[0]
    S_d_indices_broad=np.where(  (z < z_interval_d[1] +0.5) &  (z > z_interval_d[0]) )[0]
    
    S_u_indices=np.where(  (z < z_interval_u[0] ) &  (z > z_interval_u[1]) )[0]
    S_d_indices=np.where(  (z < z_interval_d[1] ) &  (z > z_interval_d[0]) )[0]
    
    S2_u_indices=np.where( (z < z_interval_u[1] ) &  (z > z_interval_u[2]) )[0] 
    S2_d_indices=np.where( (z < z_interval_d[2] ) &  (z > z_interval_d[1]) )[0]
    
    
    #H
    H=positions[H_indices]
   
    #O 
    O= positions[list(set(O_indices).intersection(set(water_indices)))]     #water
    
    OS_u= positions[list(set(O_indices).intersection(set(S_u_indices_broad)))]    #surface up
    OS_d= positions[list(set(O_indices).intersection(set(S_d_indices_broad)))]    #surface down
    #Fe
    FeS_u= positions[list(set(Fe_indices).intersection(set(S_u_indices)))]    #surface up
    FeS_d= positions[list(set(Fe_indices).intersection(set(S_d_indices)))]    #surface down
    FeS2_u= positions[list(set(Fe_indices).intersection(set(S2_u_indices)))]    #2 surface up
    FeS2_d= positions[list(set(Fe_indices).intersection(set(S2_d_indices)))]    #2 surface down
    
    #checks
    #if len(H) != N_water*2:
    #  print(f"H Err:  {len(H) - N_water*2} atom" )
    #  lost_atom_bool=True
    #if len(O) != N_water:
    #  print(f"O Err: {len(O) - N_water} atom" ) 
    #  lost_atom_bool=True
    #if len(OS_u) != N_OS:
    #  print(f"OS_u Err: {len(OS_u) - N_OS} atom" )
    #  lost_atom_bool=True
    #if len(OS_d) != N_OS:
    #  print(f"OS_d Err: {len(OS_d) - N_OS} atom" )
    #  lost_atom_bool=True
    #if len(FeS_u) != N_FeS:
    #  print(f"FeS_u Err: {len(FeS_u) - N_FeS} atom" )
    #if len(FeS_d) != N_FeS:
    #  print(f"FeS_d Err: {len(FeS_d) - N_FeS} atom" )
    #if len(FeS2_u) != N_FeS2:
    #  print(f"FeS2_u Err: {len(FeS2_u) - N_FeS2} atom" )
    #if len(FeS2_d) != N_FeS2:
    #  print(f"FeS2_d Err: {len(FeS2_d) - N_FeS2} atom" )
      
   
    return H, O, OS_u, OS_d, FeS_u, FeS_d, FeS2_u, FeS2_d, lost_atom_bool
  
def parser_low_coverage(t_slice,  N_water, N_O, N_OS, N_Fe, N_FeS, N_FeS2):
    
    H_u=[]
    H_d=[]
    O_u=[]
    O_d=[]
    OS_u=[]
    OS_d=[]
    FeS_u=[]
    FeS_d=[]
    FeS2_u=[]
    FeS2_d=[]
        
    positions=t_slice.get_positions(wrap=True)
   
    #indices species 
    H_indices = np.where(t_slice.get_atomic_numbers()==1)[0]
    O_indices = np.where(t_slice.get_atomic_numbers()==2 )[0]
    Fe_indices = np.where(t_slice.get_atomic_numbers()==3 )[0]
    
    #sort
    #print(  positions[H_indices])
    #print(  positions[H_indices_sort] [:N_water]     )
    
    #sorting
    #H
    
    H_indices_sort=np.argsort( positions[H_indices,2]) [::-1]
    H_u= positions[H_indices_sort] [:N_water]
    H_d= positions[H_indices_sort] [N_water:]
    #O 
    #water_u
    O_indices_sort=N_H+np.argsort( positions[O_indices,2]) [::-1]
    O_u=  positions[O_indices_sort][:int(N_water/2)]
    OS_u= positions[O_indices_sort][int(N_water/2):int(N_water/2)+N_OS]
    OS_d= positions[O_indices_sort][N_O -int(N_water/2) -N_OS:  N_O -int(N_water/2)   ]
    O_d=  positions[O_indices_sort][  N_O -int(N_water/2) :  ]
   
    #Fe
    Fe_indices_sort=N_H + N_O + np.argsort( positions[Fe_indices,2]) [::-1]
    FeS_u=   positions[Fe_indices_sort][:N_FeS]                       #surface up
    FeS2_u=  positions[Fe_indices_sort][N_FeS: N_FeS+ N_FeS2 ]        #2 surface up
    FeS2_d=  positions[Fe_indices_sort][ N_Fe - N_FeS- N_FeS2 : N_Fe - N_FeS ]    #2 surface down
    FeS_d=   positions[Fe_indices_sort][N_Fe - N_FeS :]   
    
  
    #no sorting (already sorted)
    #H
    #H_u=positions[H_indices[:N_water]]
    #H_d=positions[H_indices[N_water:]]
    ##O 
    ##water_u
    #O_u= positions[O_indices[:int(N_water/2)]] 
    #OS_u= positions[O_indices[int(N_water/2):int(N_water/2)+N_OS]] 
    #OS_d=positions[O_indices[ N_O -int(N_water/2) -N_OS:  N_O -int(N_water/2)   ]] 
    #O_d=positions[O_indices[  N_O -int(N_water/2) :  ]] 
    ##Fe
    #FeS_u= positions[Fe_indices[:N_FeS]]    #surface up
    #FeS2_u= positions[Fe_indices[N_FeS: N_FeS+ N_FeS2 ]]     #2 surface up
    #FeS2_d= positions[Fe_indices[ N_Fe - N_FeS- N_FeS2 : N_Fe - N_FeS ] ]   #2 surface down
    #FeS_d= positions[Fe_indices[N_Fe - N_FeS :]]
    
    #checks
    if len(H_u) != N_water:
      print(f"H Err:  {len(H_u) - N_water} atom" )
    if len(H_d) != N_water:
      print(f"H Err:  {len(H_d) - N_water} atom" )
    if len(O_u) != N_water/2:
      print(f"O Err: {len(O_u) - N_water/2} atom" ) 
    if len(O_d) != N_water/2:
      print(f"O Err: {len(O_d) - N_water/2} atom" ) 
    if len(OS_u) != N_OS:
      print(f"OS_u Err: {len(OS_u) - N_OS} atom" )
    if len(OS_d) != N_OS:
      print(f"OS_d Err: {len(OS_d) - N_OS} atom" )
    if len(FeS_u) != N_FeS:
      print(f"FeS_u Err: {len(FeS_u) - N_FeS} atom" )
    if len(FeS_d) != N_FeS:
      print(f"FeS_d Err: {len(FeS_d) - N_FeS} atom" )
    if len(FeS2_u) != N_FeS2:
      print(f"FeS2_u Err: {len(FeS2_u) - N_FeS2} atom" )
    if len(FeS2_d) != N_FeS2:
      print(f"FeS2_d Err: {len(FeS2_d) - N_FeS2} atom" )
     
    return H_u, H_d, O_u, O_d, OS_u, OS_d, FeS_u, FeS_d, FeS2_u, FeS2_d
  
def parser_water_bulk(t_slice):
    
    positions=t_slice.get_positions(wrap=True)
 
    H_indices = np.where(t_slice.get_atomic_numbers()==1)[0]
    O_indices = np.where(t_slice.get_atomic_numbers()==2 )[0]
    
    #H
    H=positions[H_indices]
   
    #O 
    O= positions[O_indices]     #water
    
   
    return H, O

# %%
#takes: atoms in a layer
#determine exactly the surfaces positions with average of z of all atom in the surface: 
def get_layer_z(A):
    z=A[:,2]
    return np.mean(z), max(z), min(z)


# %%
#FULL COVERAGE
#give #H or #O of water, z of the surface up and down , and max distance from surface
# 1)split atom above and beyond surfaces up to max distance distance. 2) (surface u) shifting z to zero, (surface d) shifting  z to zero  and inverting
def surface_up_down_split(A, z_u, z_d,  max_distance):
    buffer=0
    tmp_u=[]
    tmp_d=[]
  
    for atom in A:
        if  z_u - buffer  <= atom[2] <  z_u + max_distance  :
            tmp_u.append( [atom[0], atom[1],    atom[2]-z_u ])
            
        elif z_d - max_distance <= atom[2]<  z_d+buffer:
            tmp_d.append([atom[0], atom[1], -atom[2]+z_d])  
            
    return np.array(tmp_u), np.array(tmp_d)

            
#full_big, surface up is like the surface down but inverted with respect to all the other surfaces
def invert_raws_direction(A, L_y):
    atom_inverted=[]
  
    for atom in A:
        atom_inverted.append([atom[0], L_y- atom[1], atom[2]])  
            
    return np.array(atom_inverted)  
    
#FULL COVERAGE WITH INVERTED DOWN SURFACE (L_z=70)
#y goes to L_y-y
def surface_up_down_split_down_correction(A, z_u, z_d, L_y, max_distance):
    buffer=0
    tmp_u=[]
    tmp_d=[]
  
    for atom in A:
        if  z_u - buffer  <=atom[2] <  z_u + max_distance  :
            tmp_u.append( [atom[0], atom[1],    atom[2]-z_u ])
        elif z_d -  max_distance <= atom[2]<  z_d+buffer:
            tmp_d.append([atom[0], L_y- atom[1], -atom[2]+z_d])   #y to -y it align the oct raws
            
    return np.array(tmp_u), np.array(tmp_d)


#FOR DEFINITE NUMBER OF PARTICLES UP AND DOWN (low coverage middle coverage, or surface elements OS_u .. FeS2_d)
#does:  (surface u) shifting z to zero, (surface d) shifting  z to zero  and inverting
def surface_up_translate(A_u, z_u):
    
    tmp_u=[]
  
    for atom in A_u:
            tmp_u.append( [atom[0], atom[1],   atom[2]-z_u ])
            
    return np.array(tmp_u)

def surface_down_translate_invert(A_d, z_d ):

    tmp_d=[]
            
    for atom in A_d:    
        tmp_d.append([atom[0], atom[1], -atom[2]+z_d])
        
    return  np.array(tmp_d)

def surface_down_translate_invert_down_correction(A_d, L_y, z_d ):

    tmp_d=[]
            
    for atom in A_d:    
        tmp_d.append([atom[0],  L_y- atom[1], -atom[2]+z_d])
        
    return  np.array(tmp_d)


# %%
def write_h5py_file(all_coordinates_u, t, time_digit_lenght, h5py_files_folder):
    t_string=str(t).rjust(time_digit_lenght, ' ')
    for type, atom_positions in all_coordinates_u.items():
        with h5py.File(h5py_files_folder+type+'.h5py', 'a') as file:
            file.create_dataset(t_string, data= atom_positions, shape=atom_positions.shape, dtype=atom_positions.dtype)

# %%
#full coverage big   
def full_big_main(file_trj, L_x, L_y,  N_water, N_OS, N_FeS, N_FeS2,  z_max, z_interval_u, z_interval_d,  time_digit_lenght, h5py_files_folder):
     t=0             
     for t_slice in aio.iread(file_trj,format='lammps-dump-text'):
          #parsing
          H, O,  OS_u, OS_d, FeS_u, FeS_d, FeS2_u, FeS2_d, lost_atom_bool = parser_full_coverage(t_slice,  z_interval_u, z_interval_d,  N_water, N_OS, N_FeS, N_FeS2)

          #if lost_atom_bool:
          #     continue

          #locate the surface
          S_u, S_u_min, S_u_max =get_layer_z(OS_u)
          #S2_u, S2_u_min, S2_u_max =get_layer_z(FeS2_u)
          
          S_d, S_d_min, S_d_max =get_layer_z(OS_d) 
          #S2_d, S2_d_min, S2_d_max =get_layer_z(FeS2_d)
          #updating intervals: NOT WORKING
          #z_interval_u=[S_u_max+1,  S2_u_max +0.5 , S2_u_min - 0.8 ]
          #z_interval_d=[S_d_min-1,  S2_d_min-0.5, S2_d_max + 0.8]
         
          # translate to 0 the surfaces (and invert the d) 
          #water
          O_u, O_d=surface_up_down_split(O, S_u, S_d, z_max)   # down correction: y is inverted to align the down raws
          H_u, H_d=surface_up_down_split(H, S_u, S_d, z_max)       
          #surface atoms
          OS_u  =surface_up_translate(OS_u  , S_u)
          FeS_u =surface_up_translate(FeS_u , S_u)
          FeS2_u=surface_up_translate(FeS2_u, S_u)

          OS_d  =surface_down_translate_invert(OS_d,  S_d)
          FeS_d =surface_down_translate_invert(FeS_d, S_d)
          FeS2_d=surface_down_translate_invert(FeS2_d,S_d)

          all_coordinates={'H_u': H_u, 'O_u' : O_u,  'OS_u' : OS_u , 'FeS_u' : FeS_u , 'FeS2_u' : FeS2_u  ,
                             'H_d': H_d, 'O_d' : O_d,  'OS_d' : OS_d,  'FeS_d' : FeS_d , 'FeS2_d' : FeS2_d  }

          # wrap possible out of box in xy direction
          for key in all_coordinates.keys():
               all_coordinates[key]=invert_raws_direction(all_coordinates[key], L_y)
               wrap_boundary_xy(all_coordinates[key], L_x, L_y)
                    
          #writing h5py files 
          write_h5py_file(all_coordinates, t, time_digit_lenght, h5py_files_folder)
     
          t+=1
   
#full coverage 
def full_main(file_trj, L_x, L_y,  N_water, N_OS, N_FeS, N_FeS2,  z_max, z_interval_u, z_interval_d,  time_digit_lenght, h5py_files_folder):
     t=0             
     for t_slice in aio.iread(file_trj,format='lammps-dump-text'):
          #parsing
          H, O,  OS_u, OS_d, FeS_u, FeS_d, FeS2_u, FeS2_d , lost_atom_bool= parser_full_coverage(t_slice,  z_interval_u, z_interval_d,  N_water, N_OS, N_FeS, N_FeS2)

          if lost_atom_bool:
               continue

          #locate the surface
          S_u, S_u_min, S_u_max =get_layer_z(OS_u)
          S_d, S_d_min, S_d_max =get_layer_z(OS_d) 

          # translate to 0 the surfaces (and invert the d) 
          #water
          O_u, O_d=surface_up_down_split_down_correction(O, S_u, S_d, L_y, z_max)   # down correction: y is inverted to align the down raws
          H_u, H_d=surface_up_down_split_down_correction(H, S_u, S_d, L_y, z_max)       
          #surface atoms
          OS_u  =surface_up_translate(OS_u  , S_u)
          FeS_u =surface_up_translate(FeS_u , S_u)
          FeS2_u=surface_up_translate(FeS2_u, S_u)

          OS_d=surface_down_translate_invert_down_correction(OS_d, L_y, S_d)
          FeS_d =surface_down_translate_invert_down_correction(FeS_d, L_y,  S_d)
          FeS2_d=surface_down_translate_invert_down_correction(FeS2_d, L_y, S_d)
          #dictionary
          all_coordinates_u={'H_u': H_u, 'O_u' : O_u,  'OS_u' : OS_u , 'FeS_u' : FeS_u , 'FeS2_u' : FeS2_u  }
          all_coordinates_d={'H_d': H_d, 'O_d' : O_d,  'OS_d' : OS_d,  'FeS_d' : FeS_d , 'FeS2_d' : FeS2_d  }

          # wrap possible out of box in xy direction
          for atom_positions in all_coordinates_u.values():
               wrap_boundary_xy(atom_positions, L_x, L_y)
          for atom_positions in all_coordinates_d.values():
               wrap_boundary_xy(atom_positions, L_x, L_y)

          #writing h5py files 
          write_h5py_file(all_coordinates_u, t, time_digit_lenght, h5py_files_folder)
          write_h5py_file(all_coordinates_d, t, time_digit_lenght, h5py_files_folder)
          print(t)
          t+=1

#low coverage and middle coverage
def low_main(file_trj, L_x, L_y,  N_water, N_O, N_Fe,  N_OS, N_FeS, N_FeS2,   time_digit_lenght, h5py_files_folder):
     t=0
     for t_slice in aio.iread(file_trj,format='lammps-dump-text'):
          #parsing
          H_u, H_d, O_u, O_d, OS_u, OS_d, FeS_u, FeS_d, FeS2_u, FeS2_d = parser_low_coverage(t_slice,  N_water, N_O, N_OS, N_Fe, N_FeS, N_FeS2)

          #locate the surface
          S_u, S_u_min, S_u_max =get_layer_z(OS_u)
          S_d, S_d_min, S_d_max =get_layer_z(OS_d) 

          # translate to 0 the surfaces (and invert the d) 
          H_u   =surface_up_translate(H_u   , S_u)
          O_u   =surface_up_translate(O_u   , S_u)
          OS_u  =surface_up_translate(OS_u  , S_u)
          FeS_u =surface_up_translate(FeS_u , S_u)
          FeS2_u=surface_up_translate(FeS2_u, S_u)

          H_d   =surface_down_translate_invert_down_correction(H_d,   L_y, S_d )
          O_d   =surface_down_translate_invert_down_correction(O_d,   L_y, S_d )
          OS_d  =surface_down_translate_invert_down_correction(OS_d,  L_y, S_d)
          FeS_d =surface_down_translate_invert_down_correction(FeS_d, L_y, S_d)
          FeS2_d=surface_down_translate_invert_down_correction(FeS2_d,L_y, S_d)

          all_coordinates_u={'H_u': H_u, 'O_u' : O_u,  'OS_u' : OS_u , 'FeS_u' : FeS_u , 'FeS2_u' : FeS2_u  }
          all_coordinates_d={'H_d': H_d, 'O_d' : O_d,  'OS_d' : OS_d,  'FeS_d' : FeS_d , 'FeS2_d' : FeS2_d  }

          # wrap possible out of box in xy direction
          for atom_positions in all_coordinates_u.values():
               wrap_boundary_xy(atom_positions, L_x, L_y)
               
          for atom_positions in all_coordinates_d.values():
               wrap_boundary_xy(atom_positions, L_x, L_y)

          #writing h5py files 
          write_h5py_file(all_coordinates_u, t, time_digit_lenght, h5py_files_folder)
          write_h5py_file(all_coordinates_d, t, time_digit_lenght, h5py_files_folder)

          t+=1
                  
#down surface does not need correction!
def low_main_1_2_replicated(file_trj, L_x, L_y,  N_water, N_O, N_Fe,  N_OS, N_FeS, N_FeS2,   time_digit_lenght, h5py_files_folder):
     t=0
     for t_slice in aio.iread(file_trj,format='lammps-dump-text'):
          #parsing
          H_u, H_d, O_u, O_d, OS_u, OS_d, FeS_u, FeS_d, FeS2_u, FeS2_d = parser_low_coverage(t_slice,  N_water, N_O, N_OS, N_Fe, N_FeS, N_FeS2)

          #locate the surface
          S_u, S_u_min, S_u_max =get_layer_z(OS_u)
          S_d, S_d_min, S_d_max =get_layer_z(OS_d) 

          # translate to 0 the surfaces (and invert the d) 
          H_u   =surface_up_translate(H_u   , S_u)
          O_u   =surface_up_translate(O_u   , S_u)
          OS_u  =surface_up_translate(OS_u  , S_u)
          FeS_u =surface_up_translate(FeS_u , S_u)
          FeS2_u=surface_up_translate(FeS2_u, S_u)

          H_d   =surface_down_translate_invert(H_d,    S_d )
          O_d   =surface_down_translate_invert(O_d,    S_d )
          OS_d  =surface_down_translate_invert(OS_d,   S_d)
          FeS_d =surface_down_translate_invert(FeS_d,  S_d)
          FeS2_d=surface_down_translate_invert(FeS2_d, S_d)

          all_coordinates_u={'H_u': H_u, 'O_u' : O_u,  'OS_u' : OS_u , 'FeS_u' : FeS_u , 'FeS2_u' : FeS2_u  }
          all_coordinates_d={'H_d': H_d, 'O_d' : O_d,  'OS_d' : OS_d,  'FeS_d' : FeS_d , 'FeS2_d' : FeS2_d  }

          # wrap possible out of box in xy direction
          for atom_positions in all_coordinates_u.values():
               wrap_boundary_xy(atom_positions, L_x, L_y)
               
          for atom_positions in all_coordinates_d.values():
               wrap_boundary_xy(atom_positions, L_x, L_y)

          #writing h5py files 
          write_h5py_file(all_coordinates_u, t, time_digit_lenght, h5py_files_folder)
          write_h5py_file(all_coordinates_d, t, time_digit_lenght, h5py_files_folder)

          t+=1

# %%

#RUN IN PLACE with python3 command
path_of_run='./'

file_trj=path_of_run+"trjmagnetite.dump"
h5py_files_folder= path_of_run+"trajectory_h5py_files/" 


#ATTENITON this will delete the content of the folder
print('make_empty:',h5py_files_folder )
make_empty_folder(h5py_files_folder)


time_digit_lenght=9

f, L_x, L_y, L_z, diag_lenght, N_water, N_H, N_O, N_Fe, N_OS, N_FeS, N_FeS2= get_initial_info(file_trj)
np.save( h5py_files_folder+"diag_lenght.txt", diag_lenght)
np.save( h5py_files_folder+"scaling_factor.txt", f)

#MAIN#

if L_z >150:
     print('System == full_huge')
     #mode='full_huge'
     z_interval_u=[97.5, 95.9, 95]     
     z_interval_d=[78.4, 80 ,81] 
     z_max=78         #max distance from the surface for O huge
     np.save( h5py_files_folder+"z_max", z_max)
     full_big_main(file_trj, L_x, L_y,  N_water, N_OS, N_FeS, N_FeS2,  z_max, z_interval_u, z_interval_d,  time_digit_lenght, h5py_files_folder)
elif L_z >100:
    print('System == full_big')
    #mode='full_big'
    z_interval_u=[62.2, 60.85, 59.8]     
    z_interval_d=[43.5, 45.05, 46] 
    z_max=41          #max distance from the surface for O full_big 
    np.save( h5py_files_folder+"z_max", z_max)
    full_big_main(file_trj, L_x, L_y,  N_water, N_OS, N_FeS, N_FeS2,  z_max, z_interval_u, z_interval_d,  time_digit_lenght, h5py_files_folder)
elif L_z>70:
    print('System == full')
    #mode='full'
    z_interval_u=[42.5, 41.1, 40.1]  
    z_interval_d=[26.29, 27.45, 28.5]
    z_max=26          #max distance from the surface for O full
    np.save( h5py_files_folder+"z_max", z_max)
    full_main(file_trj, L_x, L_y,  N_water, N_OS, N_FeS, N_FeS2,  z_max, z_interval_u, z_interval_d,  time_digit_lenght, h5py_files_folder)
else:
    #mode='low'
    if f >=5:
        print('System == low, 3-6 replicated') #invert down surface
        low_main(file_trj, L_x, L_y,  N_water, N_O, N_Fe,  N_OS, N_FeS, N_FeS2,   time_digit_lenght, h5py_files_folder)
    else:
        print('System == low, 1-2 replicated') #correct down surface
        low_main_1_2_replicated(file_trj, L_x, L_y,  N_water, N_O, N_Fe,  N_OS, N_FeS, N_FeS2,   time_digit_lenght, h5py_files_folder)


# %%
