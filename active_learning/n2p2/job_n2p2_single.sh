#!/bin/bash
#SBATCH -J n2p2
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --hint=nomultithread
#SBATCH --ntasks-per-node=24
#SBATCH --qos=zen3_0512
#SBATCH --partition=zen3_0512
#SBATCH --time=600

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

#suppose you are in run_i    

echo "Working"
    mpirun -np 24 ~/n2p2-validation/bin/nnp-scaling 500 
    echo "Making dir train"
    mkdir train
    cp input.data input.nn scaling.data validation.data train
    cd train
            mpirun -np 24  ~/n2p2-validation/bin/nnp-train 
            cd ../    


~                                                                                                                                                                                                          
~                                                                                                                                                                                                          
~                                                                                                                                                                                                          
~                                                                                                                                                                                                          
~                        
