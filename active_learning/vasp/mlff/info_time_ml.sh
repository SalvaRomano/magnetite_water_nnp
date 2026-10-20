#!/bin/bash                                                                     
#RUN THIS SCRIPT IN FOLDER 38Fe/44Fe


echo "DIR     steps/fs        DFT     ML"
for k in *; do
cd $k;
for ((i=1; i <= 1000; i++)); do
DIR=run_$i
if [ -d "$DIR" ]; then
    grep E= run_$i/OSZICAR | wc -l >> step_in_run.temp 
    grep mag= run_$i/OSZICAR | wc -l >> dft_in_run.temp 
else
    break
fi
done
echo "$k"      |tr '\n' '\t'
cat step_in_run.temp | awk   '{ sum+=$1} END {print sum, sum*1.5}'  | tee -a ../step_tot.temp |tr '\n' '\t'
cat dft_in_run.temp  | awk   '{ sum+=$1} END {print sum}'           | tee -a ../dft_tot.temp       |tr '\n' '\t'
paste ../step_tot.temp ../dft_tot.temp | tail -1 | awk '{print ($1-$3)}'                         

rm step_in_run.temp
rm dft_in_run.temp

cd ../
done
cat step_tot.temp | awk  '{ sum+=$1} END {print  sum, sum*1.5, "TOTAL steps&fs"}'| tee -a sum_tot_step.tmp; rm step_tot.temp
cat dft_tot.temp | awk  '{ sum+=$1} END {print  sum, "TOTAL DFT"}'               |tee -a sum_tot_dft.tmp ; rm dft_tot.temp
paste sum_tot_step.tmp sum_tot_dft.tmp | tail -1 | awk '{print ($1-$5), "TOTAL ML"}'

rm sum_tot_step.tmp
rm sum_tot_dft.tmp


