#!/bin/bash

#SBATCH --cluster=mpp2
#SBATCH --nodes=3
#SBATCH --ntasks=64
#SBATCH --mail-type=end
#SBATCH --mail-user=stefan.haas@drei.at
#SBATCH --time=03:00:00
#SBATCH --output=Mpi3NodeTests_%A.log

source /home/hpc/t1221/lu26xag/modules
cd /home/hpc/t1221/lu26xag/cfd3d/build/

echo 64 4 4 4
mpirun -np 64 ./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver mpi --numproc 4 4 4 --linsolver jacobi --output false
