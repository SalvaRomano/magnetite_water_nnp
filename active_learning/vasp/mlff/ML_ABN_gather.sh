#!/bin/bash                                                                     
#RUN THIS SCRIPT IN FOLDER 38Fe/44Fe

for k in *; do
mkdir DATA/$k
cd $k;
for ((i=1; i <= 1000; i++)); do
DIR=run_$i
if [ -d "$DIR" ]; then
    #echo "$k/run_$i exists"
    sleep 0.001
else
    #echo "$k/run_$i does not exist"
    cp run_$(($i-1 ))/ML_ABN ../DATA/$k/
    cd ../
    break
fi
done
done
