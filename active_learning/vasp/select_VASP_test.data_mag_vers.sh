#!/bin/bash


#INPUT DATA GOOD
for f in */; do
    cd $f
    
    echo "building in $f  new test.data and append in bad_mag_water.data bad_mag_iron.data"
    if  [ -e "test_learning" ]; then
    echo "test learning present in $f"
    cd test_learning
    for k in $(for i in $(find . -mindepth 1 -type d | sort ); do cd $i; if [ -f "OUTCAR" ] && [ ! -f "bad_mag_water.txt" ] && [ ! -f "bad_mag_iron.txt" ]; then  echo $i; fi;  cd ../;  done); do
        cd $k
        if [  -f "input.data" ]; then
            echo "alles gut in $k: append in ../../test.data"    
            cat input.data >> ../../test.data
        fi    
        cd ../
    done
    for k in $(for i in $(find . -mindepth 1 -type d | sort ); do cd $i; if [ -f "bad_mag_water.txt" ]; then  echo $i; fi;  cd ../;  done); do
        cd $k
        if [  -f "input.data" ]; then
            echo "water mag in $k: append in bad_mag_water.data"    
            cat input.data >> ../../bad_mag_water.data
        fi
        cd ../
    done
    for k in $(for i in $(find . -mindepth 1 -type d | sort ); do cd $i; if [ -f "bad_mag_iron.txt" ]; then  echo $i; fi;  cd ../;  done); do
        cd $k
        if [  -f "input.data" ]; then
            echo "iron mag in $k: append in bad_mag_water.data"    
            cat input.data >> ../../bad_mag_iron.data
        fi
        cd ../
    done
    cd ../
    fi
    cd ../
done
