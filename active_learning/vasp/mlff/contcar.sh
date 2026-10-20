#!bin/bash
NH=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $1}')
NO=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $2}' )
NW=$((NH/2))
NFe=$(grep -A 1 H POSCAR | tail -1 |  awk '{print $3}' )  
#echo "$NH, $NO_TOT, $NW, $NFe" 
head -n $((8 + $NH + $NO + $NFe )) CONTCAR > temp
mv temp CONTCAR
