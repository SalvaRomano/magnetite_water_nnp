#!/bin/bash                                                                     
#SBATCH -J block5
#SBATCH -N 3                                                    
#SBATCH --partition=zen3_0512
#SBATCH --qos=zen3_0512
#SBATCH --ntasks-per-node=64                                                 
#SBATCH --ntasks-per-core=1       
#SBATCH --time=02:00:00

module purge

module load intel-oneapi-tbb/2021.7.1-gcc-12.2.0-3vlgyyf
module load intel-oneapi-compilers/2023.1.0-gcc-8.5.0-6y5tcji
module load intel-oneapi-mpi/2021.7.1-gcc-8.5.0-jduvymn
module load mkl/2022.0.1
export I_MPI_PIN_RESPECT_CPUSET=0
export OMP_NUM_THREADS=1

mpirun -np 192 /gpfs/data/fs71796/pablo/pablo/vasp_2023_04_03/bin/vasp_std

