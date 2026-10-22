#!/bin/bash
#SBATCH -N 2
#SBATCH -J LAMMPS
##SBATCH --ntasks-per-node=64
#SBATCH --partition=zen3_0512
#SBATCH --qos=zen3_0512
#SBATCH --ntasks-per-core=12
##SBATCH --time=800


module purge            # recommended

module load --auto amdscalapack/4.0-gcc-12.2.0-7neoxaq
module load --auto  amdlibflame/4.0-gcc-12.2.0-gkln5gq
module load eigen/3.4.0-gcc-12.2.0-imdheey
module load --auto gcc/12.2.0-gcc-9.5.0-ohbahza

mpirun -np 256 /home/fs71796/salvrom/n2p2-committee-modified/bin/lmp_mpi -in md.lmp

