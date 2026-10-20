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

path = os.getcwd() + "/"

max_Force_errors = []
mean_Force_errors = []

print(path)

#IMPORTANT FOR BULK WATER  add "sort = True" to "ase.io.write" (last line)
Delta_t = 50       #IMPORTANT tos et EQUAL TO THE DUMP value  in md.lmp FILE!!!!!!!  dump equi all custom 100 equi.dump id type x y z ix iy iz
also_include = []   # INCLUDE MANUALLY A CONFIG GIVEN THE TIMESTEP

u_ase = ase.io.read(path + "equi.dump", index=':', format = 'lammps-dump-text')


print("Loaded: ", path+ "equi.dump")
print(u_ase[0])

atomic_numbers = u_ase[0].get_atomic_numbers()
#oscar_order = np.argsort(atomic_numbers)
chemical_dict= {1:'H', 2:'O', 3:'Fe'}
chemical_symbols = []
for i in range(len(atomic_numbers)):
    chemical_symbols.append(chemical_dict[atomic_numbers[i]])




disagreement=np.transpose(np.loadtxt("./pruned_disagreement.txt"))

timesteps=[int(i) for i in disagreement]
print(timesteps)

# %%


try:
    os.mkdir(path + "active_learn_POSCARs/")
except:
    pass

for ts in timesteps:
    print(ts)
    try:
       os.mkdir(path + "active_learn_POSCARs/Force_Dis_config_" + str(ts).zfill(10))
    except:
       pass
    finally:
       u_ase[int(ts/Delta_t)].set_chemical_symbols(chemical_symbols)
       ase.io.write(path + "active_learn_POSCARs/Force_Dis_config_" + str(ts).zfill(10) + "/POSCAR", u_ase[int(ts/Delta_t)] )   #IMPORTANT add "sort = True"  for bulk water!


# %%
