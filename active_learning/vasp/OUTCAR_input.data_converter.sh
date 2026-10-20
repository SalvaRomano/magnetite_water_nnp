#!/bin/bash

#convert OUTCAR in input.data

for f in */; do
cd $f
rm input.data 2>/dev/null
for k in $(for i in $(find . -mindepth 2 -type d | sort )  ; do cd  $i; if [ -f "OUTCAR" ]; then  echo $i; fi;  cd ../../;  done  ); do
    cd $k
    if [ ! -f "input.data" ]; then
        echo "Working in $f/$k"
        ~/n2p2-validation/tools/converter/convert-VASP_OUTCAR.py OUTCAR input.data
        if [[ $? -ne 0 ]]; then
            echo "ERROR"
        else
            echo "Converted OUTCAR in input.data"    
        fi    
    fi
    cd ../../
    done
cd ../
done
