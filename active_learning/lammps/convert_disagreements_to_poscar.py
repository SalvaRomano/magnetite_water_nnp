# %%
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import re

import numpy as np
#import MDAnalysis
#import MDAnalysis.analysis.distances
#import MDAnalysis.lib.distances
#import MDAnalysis.transformations
#import MDAnalysis.lib.mdamath


import ase
import ase.io

# %%
#os.chdir(os.getcwd()+"/2x2x2_1300K_0/")

path = os.getcwd() + "/"

max_Force_errors = []
mean_Force_errors = []

print(path)

Delta_t = 20   #EQUAL TO THE DUMP line in  md.lmp FILE!!!!!!!  dump equi all custom 10 equi.dump id type x y z ix iy iz
also_include = []   # INCLUDE MANUALLY A CONFIG GIVEN THE TIMESTEP
t_interval=100     #minimum time interval for including structures

u_ase = ase.io.read(path + "equi.dump", index=':', format = 'lammps-dump-text')


print("Loaded: ", path+ "equi.dump")
print(u_ase[0])
#
#symbols = [128*["O"], 13*["K"], 4*["Na"], 16*["Al"], 48*["Si"]]
#symbols = [item for sublist in symbols for item in sublist]

atomic_numbers = u_ase[0].get_atomic_numbers()
#oscar_order = np.argsort(atomic_numbers)
chemical_dict= {1:'H', 2:'O', 3:'Fe'}
chemical_symbols = []
for i in range(len(atomic_numbers)):
    chemical_symbols.append(chemical_dict[atomic_numbers[i]])


#for filename in dirFiles:
#    tmp = filename.split('-')
#    ts = int(tmp[0][5:])
    #lines = open(filename, "r").readlines()
    #result = []
 #   print(ts)
    #for x in lines:
    #    result.append(float(x.split()[1]))
    #max_Force_errors.append(max(result))
    #mean_Force_errors.append(np.mean(result))
  #  timesteps.add(ts)
#lines = open("ML_ABN", "r").readlines()
#max_Force_errors = np.array(max_Force_errors)
#mean_Force_errors = np.array(mean_Force_errors)
#timesteps = np.array(timesteps)


disagreement=np.transpose(np.loadtxt("./disagreement.txt"))

timesteps_float=list(set(disagreement[0]))
timesteps=[int(i) for i in timesteps_float]
timesteps.sort()
print(timesteps)

# %%

redo = []

for i in range(len(timesteps)):
    #if (timesteps[i] - Delta_t) in redo or (timesteps[i] - 2*Delta_t) in redo or (timesteps[i] - 3*Delta_t) in redo or (timesteps[i] - 4*Delta_t) in redo or (timesteps[i] - 5*Delta_t) in redo or (timesteps[i] - 6*Delta_t) in redo:
#    gen = timesteps[i] - np.arange(1, 100000)
    gen = timesteps[i] - np.arange(1, t_interval)
#   gen = timesteps[i] - np.arange(1, 500)
#   gen = timesteps[i] - np.arange(1, 999)
    if not any(x in redo for x in gen):
        redo.append(timesteps[i])
        print(timesteps[i])

redo = redo + also_include


print(redo)
print(len(redo))


try:
    os.mkdir(path + "active_learn_POSCARs/")
except:
    pass

for ts in redo:
    print(ts)
    try:
       os.mkdir(path + "active_learn_POSCARs/Force_Dis_config_" + str(ts).zfill(10))
    except:
       pass
    finally:
       u_ase[int(ts/Delta_t)].set_chemical_symbols(chemical_symbols)
       ase.io.write(path + "active_learn_POSCARs/Force_Dis_config_" + str(ts).zfill(10) + "/POSCAR", u_ase[int(ts/Delta_t)] )   #add sort = True for bulk water!




# %%
