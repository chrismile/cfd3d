#!/bin/bash

#SBATCH --cluster=mpp2
#SBATCH --nodes=2
#SBATCH --ntasks=56
#SBATCH --mail-type=end
#SBATCH --mail-user=stefan.haas@drei.at
#SBATCH --time=01:00:00
#SBATCH --output=Mpi2NodeTests_%A.log

source /home/hpc/t1221/lu26xag/modules
cd /home/hpc/t1221/lu26xag/cfd3d/build/

echo 56 7 4 2
mpirun -np 56 ./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver mpi --numproc 7 4 2 --linsolver jacobi --output false
