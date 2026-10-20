#!/bin/bash
#
#SBATCH -J Vasp_std_check
#SBATCH -N 1
#SBATCH --qos=zen3_0512
#SBATCH --partition=zen3_0512
#SBATCH --ntasks-per-node=128
#SBATCH --ntasks-per-core=1
#SBATCH --time=20:30:00




module purge
module load intel-oneapi-tbb/2021.7.1-gcc-12.2.0-3vlgyyf
module load intel-oneapi-compilers/2023.1.0-gcc-8.5.0-6y5tcji
module load intel-oneapi-mpi/2021.7.1-gcc-8.5.0-jduvymn
module load mkl/2022.0.1

export I_MPI_PIN_RESPECT_CPUSET=0
export OMP_NUM_THREADS=1


# set -x print each command that is going to be executed and -i if there is an error exit the job
#set -x    ##optional -i


#random_number=$((RANDOM % 100 + 1))
#echo "sleeping $random"
#sleep $random_number

mpirun -np 128 /gpfs/data/fs71796/salvrom/vasp-conda_wo_conda_intel/bin/vasp_std
