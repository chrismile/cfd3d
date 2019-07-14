#!/bin/bash

#SBATCH --cluster=mpp2
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --mail-type=end
#SBATCH --mail-user=stefan.haas@drei.at
#SBATCH --time=05:00:00
#SBATCH --output=mpiDomainSplittingTests_1Node_%A.log

source /home/hpc/t1221/lu26xag/modules
cd /home/hpc/t1221/lu26xag/cfd3d/build/

echo 28 7 2 2
mpirun -np 28 ./cfd3d --scenario driven_cavity_mpiDomainSplitting --solver mpi --numproc 7 2 2 --output false
echo 27 3 3 3
mpirun -np 27 ./cfd3d --scenario driven_cavity_mpiDomainSplitting --solver mpi --numproc 3 3 3 --output false
echo 30 5 3 2
mpirun -np 30 ./cfd3d --scenario driven_cavity_mpiDomainSplitting --solver mpi --numproc 5 3 2 --output false
echo 20 5 2 2
mpirun -np 20 ./cfd3d --scenario driven_cavity_mpiDomainSplitting --solver mpi --numproc 5 2 2 --output false
