#!/bin/bash                                                                     
#RUN THIS SCRIPT IN FOLDER 38Fe/44Fe

for k in *; do
cd $k;
for ((i=1; i <= 1000; i++)); do
DIR=run_$i
if [ -d "$DIR" ]; then
    #echo "$k/run_$i exists"
    sleep 0.001
else
    #echo "$k/run_$i does not exist"
    head -n 5 run_$(($i-1 ))/ML_ABN | tail -n 1 | tee -a ../ML_LOG.txt
    cd ../
    break
fi
done
done
cat ML_LOG.txt | awk '{ sum+=$1} END {print sum}'
rm ML_LOG.txt
