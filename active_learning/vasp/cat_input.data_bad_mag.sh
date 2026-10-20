#!/bin/bash


#INPUT DATA GOOD
for f in */; do
    cd $f
    rm input.data bad* test.data  2>/dev/null
    
    echo "cat in $f new input.data bad_mag_water.data bad_mag_iron.data"
    for k in $(for i in $(find . -mindepth 2 -type d | sort ); do cd $i; if [ -f "OUTCAR" ]; then  echo $i; fi;  cd ../../;  done); do
        cd $k
        if [  -f "input.data" ]   && [ ! -f "bad_mag" ] ; then
            echo "alles gut in $f/$k: append input.data in $f"    
            cat input.data >> ../../input.data
        elif [  -f "input.data" ] && [  -f "bad_mag" ]; then
            echo "bad_mag in $f/$k: append  bad_mag.data in $f"    
            cat input.data >> ../../bad_mag.data
        fi    
        cd ../../
    done
    cd ../
done
