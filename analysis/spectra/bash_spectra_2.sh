
## Set flags
   mode=1
   nmax=50000
   z_cut=7.2 #2.00 +0.25

# launch in trj_i where chunks are to compute SSVVCF-ML-mica for each chunk


for i in chunk*/; do 
	cd $i
	cd second_layer
 	
	
	~/scripts/Yair/VSFG-magnetite/src/main.x -f pos.xyz vel.xyz -nmode ${mode}  -cell 33.584 33.584 71.47838 -dt 0.004 -max ${nmax} -lag 0.5 -rc 4 -z_cut ${z_cut} -npairs_max 9000 
	cd ../../;
done
