# start with md.lmp_T_200 in the mother folder that contains the "typ of run" like water_6_uc or random_6_uc, inside will be the T_*00/run_j
for i in T_*; do  cd $i; for j in {1..4}; do mkdir run_$j; cd run_$j;   cp  ../../md.lmp_$i ./md.lmp;  ln -s  ../../initial.data ./initial.data; cd ../;  done ; cd ../; done




