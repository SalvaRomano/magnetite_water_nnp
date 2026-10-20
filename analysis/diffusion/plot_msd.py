# %%

import numpy as np
import matplotlib.pyplot as plt
import os

# %%
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
# %%
path_of_runs="/data/home/salvatorer26/projects/ferrite_oxide/production/production_runs/replicated_1_2/"

run_folder='water_6_uc_2rep/'

T_folder='T_400/'
ns_per_step=2000


msd_folder=path_of_runs +run_folder+T_folder+"MSD.npz_files/"
#save_folder=msd_folder+'/plots'
#make_empty_folder(save_folder)

# %%
#LOAD 
msd_total={}
msd_average_para={}
msd_average_ortho={}
msd_all_para={}
msd_all_ortho={}
atom_names=['H', 'O']
surfaces=['u','d']
for surface_selection in surfaces:
    print('selected surface:', surface_selection)
    for count, atom_type in enumerate(atom_names):
        print('selected atom', atom_type)


        msd_total        [atom_type+surface_selection] = np.load(msd_folder + atom_type+ '_MSD_xy_'       + surface_selection + '.npz')
        msd_average_ortho[atom_type+surface_selection] = np.load(msd_folder + atom_type+ '_MSD_xy_ortho_' + surface_selection + '.npz')
        msd_average_para [atom_type+surface_selection] = np.load(msd_folder + atom_type+ '_MSD_xy_para_'  + surface_selection + '.npz')
        msd_all_para     [atom_type+surface_selection] = np.load(msd_folder + atom_type +'_MSD_xy_all_atoms_para_'  + surface_selection +  '.npz')
        msd_all_ortho    [atom_type+surface_selection] = np.load(msd_folder + atom_type +'_MSD_xy_all_atoms_ortho_' + surface_selection +  '.npz')
        
nanosec=np.arange(0, len(msd_total['Hu']["msd"]) )/ns_per_step #same axes for every plot

# %%
#PLOT MSD

for atom in atom_names:
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))

    for i, ax in enumerate(axes):
        if i == 0:
            ax.plot(nanosec, msd_total        [atom+surfaces[i]]["msd"], label = 'tot')
            ax.plot(nanosec, msd_average_ortho[atom+surfaces[i]]["msd"], label = 'ortho')
            ax.plot(nanosec, msd_average_para [atom+surfaces[i]]["msd"], label = 'para')
            ax.set_title(atom+" MSD "+surfaces[i] )
            ax.set_xlabel("ns")
            ax.set_ylabel("$\\AA^2$")
            ax.set_xlim(0)
            ymin_1, ymax_1 = ax.set_ylim(0)
            ax.legend()
        else:             
            ax.plot(nanosec, msd_total        [atom+surfaces[i]]["msd"], label = 'tot')
            ax.plot(nanosec, msd_average_ortho[atom+surfaces[i]]["msd"], label = 'ortho')
            ax.plot(nanosec, msd_average_para [atom+surfaces[i]]["msd"], label = 'para')
            ax.set_title(atom+" MSD "+surfaces[i] )
            ax.set_xlabel("ns")
            ax.set_ylabel("$\\AA^2$")
            ax.set_xlim(0)
            ymin_2, ymax_2 = ax.set_ylim(0)
            ax.legend()

    common_ylim = (0, np.max([ymax_1,ymax_2]))
    for ax in axes:
        ax.set_ylim(common_ylim)
    plt.subplots_adjust(wspace=0.4)       
    plt.tight_layout()
    plt.show()    
    #plt.savefig(f"{save_folder}/{atom}_MSD", bbox_inches="tight", facecolor="white", transparent=True, dpi=200) 
# %%    
#PLOT PARA MSD ALL ATOM
 
for atom in atom_names:
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))


    for i, ax in enumerate(axes):
        if i == 0:
            ax.plot(nanosec, msd_all_para [atom+surfaces[i]]["msd"])
            ax.set_title(atom+" para all atoms MSD "+surfaces[i] )
            ax.set_xlabel("ns")
            ax.set_ylabel("$\\AA^2$")
            ax.set_xlim(0)
            ymin_1, ymax_1 = ax.set_ylim(0)
        else:      
            ax.plot(nanosec, msd_all_para [atom+surfaces[i]]["msd"])
            ax.set_title(atom+" para all atoms MSD "+surfaces[i] )
            ax.set_xlabel("ns")
            ax.set_ylabel("$\\AA^2$")
            ax.set_xlim(0)
            ymin_2, ymax_2 = ax.set_ylim(0)


    common_ylim = (0, np.max([ymax_1,ymax_2]))
    for ax in axes:
        ax.set_ylim(common_ylim)
    plt.subplots_adjust(wspace=0.4)       
    # Display the plot
    plt.tight_layout()
    plt.show()             
    #plt.savefig(f"{save_folder}/{atom}_MSD_all_para", bbox_inches="tight", facecolor="white", transparent=True, dpi=200)
# %%    
#PLOT ORTHO MSD ALL ATOM
 
for atom in atom_names:
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4) )


    for i, ax in enumerate(axes):
        if i == 0:
            ax.plot(nanosec, msd_all_ortho [atom+surfaces[i]]["msd"])
            ax.set_title(atom+" ortho all atoms MSD "+surfaces[i] )
            ax.set_xlabel("ns")
            ax.set_ylabel("$\\AA^2$")
            ax.set_xlim(0)
            ymin_1, ymax_1 = ax.set_ylim(0)
        else:      
            ax.plot(nanosec, msd_all_ortho [atom+surfaces[i]]["msd"])
            ax.set_title(atom+ " ortho all atoms MSD "+surfaces[i] )
            ax.set_xlabel("ns")
            ax.set_ylabel("$\\AA^2$")
            ax.set_xlim(0)
            ymin_2, ymax_2 = ax.set_ylim(0)

    common_ylim = (0, np.max([ymax_1, ymax_2]))
    for ax in axes:
        ax.set_ylim(common_ylim)
    plt.subplots_adjust(wspace=0.4)       
    # Display the plot
    plt.tight_layout()
    plt.show()             
    #plt.savefig(f"{save_folder}/{atom}_MSD_all_ortho", bbox_inches="tight", facecolor="white", transparent=True, dpi=200)

#%%
#SINGLE PLOTS
#
#
#
#
#plt.plot(nanosec, msd_total["msd"], label = 'tot')
#plt.plot(nanosec, msd_average_ortho["msd"], label = 'ortho')
#plt.plot(nanosec, msd_average_para["msd"], label = 'para')
#plt.title("MSD")
#plt.xlabel("ns")
#plt.ylabel("$\\AA^2$")
#plt.legend()
#
#save_path="/data/home/salvatorer26/Desktop/group_meeting/Pres_9_cluster"
##plt.savefig(f"{save_path}/{run_type}MSD_{surface_selection}", bbox_inches="tight", facecolor="white", transparent=True, dpi=200)
#
#plt.show()
#
##%%
#msd_all_para = np.load(msd_folder + 'O_MSD_xy_all_atoms_para'   +"_" + surface_selection +   '.npz')
#plt.plot(nanosec, msd_all_para["msd"])
#plt.title("para MSD")
#plt.xlabel("ns")
#plt.ylabel("$\\AA^2$")
#
##plt.savefig(f"{save_path}/{run_type}MSD_para_{surface_selection}", bbox_inches="tight", facecolor="white", transparent=True, dpi=200)
#
#plt.show()
##%%
#
#
#
#msd_all_ortho= np.load(msd_folder + 'O_MSD_xy_all_atoms_ortho_' + surface_selection +   '.npz')
#plt.plot(nanosec, msd_all_ortho["msd"])
#plt.title("Ortho MSD")
#plt.xlabel("ns")
#plt.ylabel("$\\AA^2$")
#
#
#
##plt.savefig(f"{save_path}/{run_type} MSD_ortho_{surface_selection}", bbox_inches="tight", facecolor="white", transparent=True, dpi=200)
#
#plt.show()
#


