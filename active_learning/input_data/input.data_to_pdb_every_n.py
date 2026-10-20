# %%
import numpy as np
from matplotlib import pyplot as plt
import re
import os

import MDAnalysis as mda
import MDAnalysis.analysis.distances
import MDAnalysis.lib.distances
import MDAnalysis.transformations
import MDAnalysis.lib.mdamath as mdamath

Qm_Calculation_Counter = 0
Qm_Calculation_list = []

lines=open("input.data", "r").readlines()

every_n=10


#os.remove("input.data_1")
frame = 0
file = open(str(frame) + '_' "input.data", 'a')
#
#for i in range(len(lines)):
#    test = lines[i].strip()
#    if test == 'begin':
#        frame += 1
#        file.write('ITEM: TIMESTEP \n')
#        file.write(str(frame) + '\n')
#        file.write('ITEM: BOX BOUNDS xy xz yz pp pp pp \n')
#    if test[0] == 'lattice':

# %%
lines = open("input.data", "r").readlines()

#os.remove("input.data_1")
#file = open("input.data", 'a')

file = open("XDATCAR.data", 'a')
frame = 0
struct_count=0
for i in range(1, len(lines)):
    test = lines[i].strip()
#    print(test[0:19])
    if test == 'begin':
        struct_count+=1
        if  struct_count % every_n == 0:
        #print("yo")
            comment = lines[i + 1].strip()
            xlo, xhi, xy = map(float, lines[i + 2].split()[1:4])
            ylo, yhi, xz = map(float, lines[i + 3].split()[1:4])
            zlo, zhi, yz = map(float, lines[i + 4].split()[1:4])

            positions = []
            elements = []
            k = 5
            #print(lines[i + k].split()[0])
            while lines[i + k].split()[0] == 'atom':
                positions.append([float(p) for p in lines[i + k].split()[1:4]])
                elements.append(lines[i + k].split()[4])
                k += 1
        #print(positions)

            u = mda.Universe.empty(len(elements), 1, trajectory= True)
            u.add_TopologyAttr('type', elements)
            u.add_TopologyAttr('name', elements)

#        u.load_new(positions, order='fac')
            u.atoms.positions = positions
        
            box = np.zeros((3, 3), dtype=np.float64)
            box[0] = xlo, xhi, xy
            box[1] = ylo, yhi, xz
            box[2] = zlo, zhi, yz

            xlen, ylen, zlen, alpha, beta, gamma = mdamath.triclinic_box(*box)
            u.dimensions = mdamath.triclinic_box(*box)

            with MDAnalysis.Writer("input_" + str(frame) + ".pdb", multiframe=False) as pdb:
                pdb.write(u)
            frame += every_n
            print(frame)
#        N_of_elem = 0
#        N_atoms = 0
#        chemistry = []
#        while lines[i + 2 + N_of_elem].strip()[0].isalpha():
#            #    print(lines[i + 2 + N_of_elem].strip())
#            element = " ".join(re.findall("[a-zA-Z]+", lines[i + 2 + N_of_elem].strip()))
#            N_per_element = int(" ".join(re.findall("[0-9]+", lines[i + 2 + N_of_elem].strip())))
#        #    print(element, N_per_element)
#            chemistry.append([element, N_per_element])
#            N_of_elem += 1
#        #print(chemistry)
#        for k in range(len(chemistry)):
#            N_atoms += chemistry[k][1]
#        Primitive_lattice_start = i + N_of_elem + 5
#        if i == 1918:
#            print(Primitive_lattice_start)
#        Atomic_positions_start = i + N_of_elem + 11
#        Total_Energy_start = i + 12 + N_atoms + N_of_elem + 2
#        Forces_start = i + 12 + N_atoms + N_of_elem + 6
#        file.write("begin" + '\n')
#        file.write("comment" + '\n')
#        for k in range(3):
#            file.write("lattice " + lines[Primitive_lattice_start + k].strip() + '\n')
#        for k in range(len(chemistry)):
#            for l in range(chemistry[k][1]):
#                file.write("atom " + lines[Atomic_positions_start + l].strip() + ' ' + chemistry[k][0] + ' '  '<c1> ' + '<n1> ' + lines[Forces_start + l].strip() + ' ' + '\n')
#            Atomic_positions_start += chemistry[k][1]  # To read in new elementtype
#            Forces_start += chemistry[k][1]  # To read in new elementtype
#
#        file.write("energy " + lines[Total_Energy_start].strip() + '\n')
#        file.write("charge " + '<unused>' + '\n')
#
#        file.write("end" + '\n')

# %%
