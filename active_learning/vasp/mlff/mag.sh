
#small = 38 Fe 56 O     MAGMOM = (2N_water)*0 (N_water+ N_oxygen)*0  4*4 3*-4 2*4 2*-4 4*4 2*-4 4*4 2*-4 4*4 2*-4 2*4 3*-4 4*4              Final Mag (M in the script) 48

#big   = 44 Fe 64 O     MAGMOM = (2N_water)*0 (N_water+ N_oxygen)*0  4*4 3*-4 2*4 2*-4 4*4 2*-4 4*4 2*-4 4*4 2*-4 4*4 2*-4 2*4 3*-4 4*4     Final MAG  (L in the script) 56


#INCAR MAGMOM ====  INITIAL MAGNETIZATION 
#         H    O     Fe
#MAGMOM = 28*0 70*0 4*4 3*-4 2*4 2*-4 4*4 2*-4 4*4 2*-4 4*4  2*-4 2*4 3*-4 4*4 



#!bin/bash
NH=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $1}')
NO_TOT=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $2}' )
NW=$((NH/2))
NFe=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $3}' )  
#echo "$NH, $NO_TOT, $NW, $NFe" 

if [[ "$NFe" -eq 38 ]]; then
    echo "$NFe"
    sed -i "s/MAGMOM/MAGMOM = $NH*0 $NO_TOT*0 4*4 3*-4 2*4 2*-4 4*4 2*-4 4*4 2*-4 4*4 2*-4 2*4 3*-4 4*4/g" INCAR
elif [[ "$NFe" == 44 ]]; then
    echo "$NFe"
    sed  -i "s/MAGMOM/MAGMOM = $NH*0 $NO_TOT*0 4*4 3*-4 2*4 2*-4 4*4 2*-4 4*4 2*-4 4*4 2*-4 4*4 2*-4 2*4 3*-4 4*4/g" INCAR
else
    echo "UNEXPECTED NUMBER OF IRON"
fi

