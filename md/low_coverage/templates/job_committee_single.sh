#!/bin/bash

#SBATCH -N 1
##SBATCH -n 12
#SBATCH -J committee
#SBATCH --partition=zen3_0512
#SBATCH --qos=zen3_0512
#SBATCH --ntasks-per-core=1
# module load <modules>
##SBATCH --time=500

module purge            # recommended

module load --auto amdscalapack/4.0-gcc-12.2.0-7neoxaq
module load --auto  amdlibflame/4.0-gcc-12.2.0-gkln5gq
module load eigen/3.4.0-gcc-12.2.0-imdheey
module load --auto gcc/12.2.0-gcc-9.5.0-ohbahza





mpirun -np 128 /home/fs71576/salvatore/n2p2-committee-modified/bin/lmp_mpi -in md.lmp 

