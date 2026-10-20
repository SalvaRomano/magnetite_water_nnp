# %%
import numpy as np
import matplotlib.pyplot as plt
import re
import h5py
import ase.io as aio



# %%
import os

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
#extract z for u and d and put in the same list
def z_all_t_extractor(H):
    zH=[]
    for atom in H:
        zH.append(atom[2])
    return zH
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
#oct and tet analysis based on oxygens surface layer: OS or OS_d as input, box xy diagonal, and  cluster parm epsilon =1

#output diagonal O raws displacment to separate (Fe) oct and tet zones
def get_diagonal_O_raws(OS, control, diag_lenght,  epsilon):
    diagonal_projection=[]

    
    for atom in OS:
        r=np.sqrt(atom[0]**2+atom[1]**2)
        theta=np.arctan(atom[1]/atom[0])
        #print(theta)
        if control == 'u':
            diagonal_projection.append(r*np.cos(theta - (np.pi/4)  )   )
        elif control == 'd':
            diagonal_projection.append(r*np.sin(theta - (np.pi/4)  )   )
        else:
            print("control must be 'u' or 'd'")
            return

    #clusterization of diagonal projection of O position in raws
    displacement=[]
    
    displacement.append([diagonal_projection[0]])

    for i, proj in enumerate(diagonal_projection[1:]):
        copy=displacement
     
        for j in range(len(copy)):
            #print(j)
            if np.abs(proj-displacement[j][0])<epsilon:
                displacement[j].append(proj)
                #print("existing")
                break
            elif j==len(copy)-1:      #last round and no match found create new entry
                #print("new")
                displacement.append([proj])


    #average over the O diagonal raws
    diag_raws=[]   
    for i in displacement:
        diag_raws.append(np.mean(i))
        
    diag_raws=list(np.sort(diag_raws))
   
    if control == 'u':
        diag_raws.insert(0,0)
        diag_raws.append(diag_lenght)
    else:
        diag_raws.insert(0,-diag_lenght/2)
        diag_raws.append(diag_lenght/2)
    #print(diag_raws)
    return diag_raws


#oct/tet regions atomic separation; inside O raws = oct, outside = tet: 
def tet_oct_split(O, diag_raw , control ):
      
    O_oct=[]
    O_tet=[]
   
    for atom in O:
        r=np.sqrt(atom[0]**2+atom[1]**2)
        
        theta=np.arctan(atom[1]/atom[0])
    
        if control == 'u':
            diagonal_projection=(r*np.cos(theta - (np.pi/4)  ))     #up cos raws   \\\\\           
        elif control == 'd':
            diagonal_projection=(r*np.sin(theta - (np.pi/4)  ))     #down sin raws  /////   
            
        else:  
            print("control must be 'u' or 'd'")
            return   
        
        if diagonal_projection <diag_raw[0] or diagonal_projection > diag_raw[-1]:  #inbound check
                print("OUT of bounds: diag =" , diagonal_projection, "cartesian =", atom)
                
        if (diag_raw[2]-diag_raw[1])<(diag_raw[3]-diag_raw[2]):           #check if the first raw (the zero is the border)is oct or tet (if T == oct if F ==tet )
            # first raw is oct
            #if diagonal projection in (diag_raw[2n],diag_raw[2n+1])   then tet region  even!
            #if diagonal projection in (diag_raw[2n+1],diag_raw[2n+2]) then oct region odd!
            for n in range(len(diag_raw)-1):  
                if diag_raw[n] <= diagonal_projection <= diag_raw[n+1]: 
                    if n % 2 ==0:
                        O_tet.append([atom[0], atom[1], atom[2]]) #n even ==> outside O raws == tet region
                    else:
                        O_oct.append([atom[0], atom[1], atom[2]]) #n odd ==> inside O raws == oct region
                    break
        else:
            #first raw is tet
            for n in range(len(diag_raw)-1):  
                if diag_raw[n] <= diagonal_projection <= diag_raw[n+1]: 
                    if n % 2 ==0:
                        O_oct.append([atom[0], atom[1], atom[2]]) #n even ==> inside O raws == oct region
                    else:
                        O_tet.append([atom[0], atom[1], atom[2]]) #n odd ==> outside O raws == tet region
                    break
        
    if len(O_oct)+len(O_tet) != len(O):
        print("MISSING TET/OCT ATOM: ",len(O)-  (len(O_oct)+len(O_tet))) 
    return O_oct,  O_tet

# %%

def denisity_hist_main(path_h5py, diag_lenght, max_dist, n_bins, surface_selection,  control ):

     #init hist
     bins = np.linspace(0, max_dist, n_bins) 
     bins_centers=(bins[1:]+bins[:-1])/2


     H_hist=[0]*(n_bins-1)
     O_hist=[0]*(n_bins-1)
     N_H=0       #Total Number of H up and down in every t
     N_O=0
     H_tet_hist=[0]*(n_bins-1)
     H_oct_hist=[0]*(n_bins-1)
     O_tet_hist=[0]*(n_bins-1)
     O_oct_hist=[0]*(n_bins-1)
     N_H_oct=0
     N_H_tet=0
     N_O_oct=0
     N_O_tet=0
     
     #PROFILE DENSITY LOOP
     #time
     t=0
     h5py_str='.h5py'
     with h5py.File(path_h5py+"H_"+surface_selection+h5py_str, 'r') as f_H, h5py.File(path_h5py+"O_"+surface_selection+h5py_str, 'r') as f_O, h5py.File(path_h5py+"OS_"+surface_selection+h5py_str, 'r') as f_OS:
         for k in f_H.keys(): #main time loop

               H = f_H[k][:]
               O = f_O[k][:]
               OS = f_OS[k][:]

               if t %10 == 0:   #it does every 10 steps
                    #get O raws info
                    diag_raw=get_diagonal_O_raws(OS, control, diag_lenght, 1 ) #give OS or OS_d, control u or d accordingly, diag_lenght and cluster param

               #split oct/tet
               O_oct, O_tet=tet_oct_split(O, diag_raw , control )     
               H_oct, H_tet=tet_oct_split(H, diag_raw , control )

               #extract z and unify
               zH=z_all_t_extractor(H)                 
               zO=z_all_t_extractor(O)
               zH_tet= z_all_t_extractor(H_tet)
               zH_oct= z_all_t_extractor(H_oct)
               zO_tet= z_all_t_extractor(O_tet)
               zO_oct= z_all_t_extractor(O_oct)

               #make instograms NOT normalized, NOT density
               H_hist+=np.histogram(zH, bins=bins)[0]      
               O_hist+=np.histogram(zO, bins=bins)[0]
               H_tet_hist+=np.histogram(zH_tet, bins=bins)[0]
               H_oct_hist+=np.histogram(zH_oct, bins=bins)[0]
               O_tet_hist+=np.histogram(zO_tet, bins=bins)[0] 
               O_oct_hist+=np.histogram(zO_oct, bins=bins)[0]   

               #Normalization: tracking the number in each time
               N_H+=len(zH)        
               N_O+=len(zO)        
               N_H_tet+=len(zH_tet)
               N_H_oct+=len(zH_oct)
               N_O_tet+=len(zO_tet)
               N_O_oct+=len(zO_oct)
           
               t+=1
               
     histograms = {
                    'H_hist': H_hist,
                    'O_hist': O_hist,
                    'H_oct_hist': H_oct_hist,
                    'O_oct_hist': O_oct_hist,
                    'H_tet_hist': H_tet_hist,
                    'O_tet_hist': O_tet_hist
               }
     
     counts = {
              'N_H': N_H,
              'N_O': N_O,
              'N_H_oct': N_H_oct,
              'N_O_oct': N_O_oct,
              'N_H_tet': N_H_tet,
              'N_O_tet': N_O_tet
          }
     

              
     return histograms, counts, bins_centers, t
               
     

# %%
def save_histograms_counts(path_hist_density, surface, histograms, counts, bins_centers, t):
    types = ['H', 'O']
    regions = ['', '_oct', '_tet']

    for type in types:
        for reg in regions:
            hist_key = type + reg + '_hist'
            count_key = 'N_' + type + reg
            filename = path_hist_density + type + reg + '_hist_' + surface

            np.save(filename, (bins_centers, histograms[hist_key], counts[count_key]))

    np.save(path_hist_density + 'tot_time', t)

# %%
#single run (HIST FOLDER "path_hist_density" must be empty or not existing)

path_of_run='./'
h5py_files_folder="trajectory_h5py_files/" 
name_hist_density_folder="/histogram_density/"


control='u'   #do not change! (because h5py are already down corrected)


path_h5py=path_of_run+h5py_files_folder
path_hist_density=path_of_run+name_hist_density_folder

make_folder(path_hist_density)

if os.path.exists(path_h5py+"diag_lenght.txt.npy"):
    diag_lenght=np.load(path_h5py+"diag_lenght.txt.npy")
if os.path.exists(path_h5py+"z_max.npy"):
    max_dist=26   #manual
    n_bins=500

print("max_dist = ",max_dist)
print("diag_lenght = ",diag_lenght)
print("n_bins =", n_bins )
    
for surface_selection in ['d','u']:      
    print("Surface selected =" , surface_selection)
    histograms, counts, bins_centers, t= denisity_hist_main(path_h5py, diag_lenght, max_dist, n_bins, surface_selection,  control )
    #Saving hist in folder
    save_histograms_counts(path_hist_density, surface_selection, histograms, counts, bins_centers, t)

# %%
#run all 
#
#path_of_runs='/data/home/salvatorer26/projects/ferrite_oxide/production/production_runs/'
#h5py_files_folder="trajectory_h5py_files" 
#name_hist_density_folder="/histogram_density/"
#control='u'   #not needed d because h5py are already down corrected
#surface_selection=['u','d']       
#
#path_h5pys=find_folder_by_name(path_of_runs, h5py_files_folder)
#
#for path_h5py in path_h5pys:
#    path_hist_density = os.path.dirname(path_h5py) + name_hist_density_folder 
#    path_h5py=path_h5py+'/'
#    
#    print(path_h5py)
#    print(path_hist_density)
#    make_empty_folder(path_hist_density)
#    
#    
#    if os.path.exists(path_h5py+"diag_lenght.txt.npy"):
#        diag_lenght=np.load(path_h5py+"diag_lenght.txt.npy")
#    if os.path.exists(path_h5py+"z_max.npy"):
#        max_dist=np.load(path_h5py+"z_max.npy")
#        if max_dist > 30:           #max dist are 8, 26 and 41
#            n_bins=int(500*41/26)
#        else:
#            n_bins=500    
#    else:
#        max_dist=8   #low coverage distance
#        n_bins=300
#    
#    print("max_dist = ",max_dist)
#    print("diag_lenght = ",diag_lenght)
#    print("n_bins =", n_bins )
#    
#    for surface in surface_selection:
#        print("Surface selected =" , surface)
#        histograms, counts, bins_centers, t= denisity_hist_main(path_h5py, diag_lenght, max_dist, n_bins, surface,  control )
#        #Saving hist in folder
#        save_histograms_counts(path_hist_density, surface, histograms, counts, bins_centers, t)

# %%

#SURFACE DENSITY LOOP

# control='u'  #no need of 'd' h5py down files are already flipped
# #control_d='d'    #second paper control_d='d'   , here is not needed because the 2 surfaces have the same orientation oct/tet 

# max_dist=41    #max cut distance from each surface: <26  for full

# #hinstograms init
# n_bin=500
# dx=max_dist/(n_bin-1)
# bins = np.linspace(0, max_dist, n_bin) 
# bins_centers=(bins[1:]+bins[:-1])/2


# H_hist=[0]*(n_bin-1)
# O_hist=[0]*(n_bin-1)
# N_H=0       #Total Number of H up and down in every t
# N_O=0
# H_tet_hist=[0]*(n_bin-1)
# H_oct_hist=[0]*(n_bin-1)
# O_tet_hist=[0]*(n_bin-1)
# O_oct_hist=[0]*(n_bin-1)
# N_H_oct=0
# N_H_tet=0
# N_O_oct=0
# N_O_tet=0

# #time
# t=0

# surface_selection='_d'
# h5py_str='.h5py'
# with h5py.File(h5py_files_folder+"H"+surface_selection+h5py_str, 'r') as f_H, h5py.File(h5py_files_folder+"O"+surface_selection+h5py_str, 'r') as f_O, h5py.File(h5py_files_folder+"OS"+surface_selection+h5py_str, 'r') as f_OS:
#     for k in f_H.keys(): #main time loop
      
#           H = f_H[k][:]
#           O = f_O[k][:]
#           OS = f_OS[k][:]
     
#           if t %10 == 0:   #it does every 10 steps
#                #get O raws info
#                diag_raw=get_diagonal_O_raws(OS, control, diag_lenght, 1 ) #give OS or OS_d, control u or d accordingly, diag_lenght and cluster param
          
#           #split oct/tet
#           O_oct, O_tet=tet_oct_split(O, diag_raw , control )     
#           H_oct, H_tet=tet_oct_split(H, diag_raw , control )
          
          
#           #extract z and unify
#           zH=z_all_t_extractor(H)                 
#           zO=z_all_t_extractor(O)
#           zH_tet= z_all_t_extractor(H_tet)
#           zH_oct= z_all_t_extractor(H_oct)
#           zO_tet= z_all_t_extractor(O_tet)
#           zO_oct= z_all_t_extractor(O_oct)

#           #make instograms NOT normalized, NOT density
#           H_hist+=np.histogram(zH, bins=bins)[0]      
#           O_hist+=np.histogram(zO, bins=bins)[0]
#           H_tet_hist+=np.histogram(zH_tet, bins=bins)[0]
#           H_oct_hist+=np.histogram(zH_oct, bins=bins)[0]
#           O_tet_hist+=np.histogram(zO_tet, bins=bins)[0] 
#           O_oct_hist+=np.histogram(zO_oct, bins=bins)[0]   
          
#           #Normalization: tracking the number in each time
#           N_H+=len(zH)        
#           N_O+=len(zO)        
#           N_H_tet+=len(zH_tet)
#           N_H_oct+=len(zH_oct)
#           N_O_tet+=len(zO_tet)
#           N_O_oct+=len(zO_oct)
#           #updates counts
#           #print(t)
#           t+=1
     

# # %%
# #save histograms
# def save_histograms(folder_hist, surface, H_hist, O_hist,  H_tet_hist, H_oct_hist, O_tet_hist, O_oct_hist, N_H, N_O, N_H_tet, N_H_oct, N_O_tet, N_O_oct, bins_centers, t):

#     types=['H', 'O']
#     hist_str='_hist'
#     regions=['', '_oct', '_tet' ]
#     for type in types:
#         for reg in regions:
#             print(type+reg+hist_str+ surface)
#                                 #name                                   #bin center         #hist               #normalization          
#             np.save(folder_hist+type+reg+hist_str+surface, ( bins_centers, eval(type+reg+hist_str)   ,  eval("N_" +type +reg)  ))
#     np.save(folder_hist+'tot_time', t)

# # %%
# #save histograms

# folder_hist="./histograms_density/"
# make_folder(folder_hist)

# types=['H', 'O']
# hist_str='_hist'
# regions=['', '_oct', '_tet' ]
# str_bins_centers='_bins_centers'
# for type in types:
#     for reg in regions:
#         print(type+reg+hist_str+ surface_selection)
#                             #name                                   #bin center         #hist               #normalization          
#         np.save(folder_hist+type+reg+hist_str+surface_selection, ( bins_centers, eval(type+reg+hist_str)   ,  eval("N_" +type +reg)  ))
# np.save(folder_hist+'tot_time', t)


