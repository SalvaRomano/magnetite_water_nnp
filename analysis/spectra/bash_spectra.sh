
## Set flags
   mode=1
   nmax=50000
   z_cut=3.70 #2.00 +0.25

# compute SSVVCF-ML-mica


for i in chunk*/; do 
	cd $i
 	~/scripts/VSFG-magnetite-main/src/main.x -f pos.xyz vel.xyz -nmode ${mode}  -cell 33.584 33.584 71.47838  -dt 0.004 -max ${nmax} -lag 0.5 -rc 4 -z_cut ${z_cut}  
	cd ../;
done
