import re
import os

Qm_Calculation_Counter = 0
Qm_Calculation_list = []

lines=open("ML_ABN", "r").readlines()

#os.remove("input.data_1")
file = open("input.data", 'a')


for i in range(1,len(lines)):
    test = lines[i].strip()
#    print(test[0:19])
    if test[0:27] == 'Atom types and atom numbers':
        N_of_elem = 0
        N_atoms = 0
        chemistry = [] 
        while lines[i + 2 + N_of_elem].strip()[0].isalpha():
        #    print(lines[i + 2 + N_of_elem].strip())
            element = " ".join(re.findall("[a-zA-Z]+", lines[i + 2 + N_of_elem].strip()))
            N_per_element = int(" ".join(re.findall("[0-9]+", lines[i + 2 + N_of_elem].strip())))
        #    print(element, N_per_element)
            chemistry.append([element, N_per_element])
            N_of_elem += 1
        #print(chemistry)
        for k in range(len(chemistry)):
            N_atoms += chemistry[k][1] 
        Primitive_lattice_start = i + N_of_elem + 5 + 4
        if i == 1918:
            print(Primitive_lattice_start)
        Atomic_positions_start = i + N_of_elem + 11 + 4
        Total_Energy_start = i + 12 + N_atoms + N_of_elem + 2 + 4
        Forces_start = i + 12 + N_atoms + N_of_elem + 6 + 4
        file.write("begin" + '\n')
        file.write("comment" + '\n')
        for k in range(3):
            file.write("lattice " + lines[Primitive_lattice_start + k].strip() + '\n')
        for k in range(len(chemistry)):
            for l in range(chemistry[k][1]):
                file.write("atom " + lines[Atomic_positions_start + l].strip() + ' ' + chemistry[k][0] + ' '  '<c1> ' + '<n1> ' +  lines[Forces_start + l].strip() + ' ' + '\n')
            Atomic_positions_start += chemistry[k][1]    #To read in new elementtype
            Forces_start += chemistry[k][1]    #To read in new elementtype

        file.write("energy " + lines[Total_Energy_start].strip() + '\n')
        file.write("charge " + '<unused>' + '\n')

        file.write("end" + '\n')
