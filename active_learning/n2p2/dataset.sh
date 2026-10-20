#!/bin/bash



module purge

module load gcc/10.2.0-gcc-9.1.0-2aa5hfe
module load ncurses/6.3-gcc-12.2.0-eg7rbxv
module load libiconv/1.16-gcc-12.2.0-matm3cb
module load python/3.9.15-gcc-12.2.0-my6jxu2
module load amdblis/4.0-gcc-12.2.0-efnclyg
module load amdlibflame/4.0-gcc-12.2.0-gkln5gq
module load openmpi/4.1.4-gcc-12.2.0-sugs3ze
module load amdscalapack/4.0-gcc-12.2.0-rfeducz
module load eigen/3.4.0-gcc-12.2.0-imdheey



export LD_LIBRARY_PATH=$LIBRARY_PATH:$LD_LIBRARY_PATH

weights=( 26 15 31 21 31  )   #change weights according to  "learn_weights.ipynb", and if it is a different system change the  number for the athomic species   "001=H" should change

index=0

for k in $(for i in $(find . -mindepth 2 -maxdepth 2 -type d | sort )  ; do echo $i; done); do
echo "$k"
cd $k
        mkdir dataset 2>/dev/null
        cp weights.001.0000${weights[$index]}.out dataset/weights.001.data
        cp weights.008.0000${weights[$index]}.out dataset/weights.008.data
        cp weights.026.0000${weights[$index]}.out dataset/weights.026.data
        cp input.data input.nn scaling.data dataset/
        cd dataset
             mpirun -np 64  ~/n2p2-validation/bin/nnp-dataset 1
        cd ../../../
        ((index++))
done
wait

