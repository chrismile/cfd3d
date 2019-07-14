#!/bin/bash

#SBATCH --cluster=mpp2
#SBATCH --nodes=2
#SBATCH --ntasks=56
#SBATCH --mail-type=end
#SBATCH --mail-user=stefan.haas@drei.at
#SBATCH --time=08:00:00
#SBATCH --output=mpiDomainSplittingTests_2Node_%A.log

source /home/hpc/t1221/lu26xag/modules
cd /home/hpc/t1221/lu26xag/cfd3d/build/

echo 56 7 4 2
mpirun -np 56 ./cfd3d --scenario driven_cavity_mpiDomainSplitting --solver mpi --numproc 7 4 2 --output false
echo 64 4 4 4
mpirun -np 64 ./cfd3d --scenario driven_cavity_mpiDomainSplitting --solver mpi --numproc 4 4 4 --output false
echo 60 5 4 3
mpirun -np 60 ./cfd3d --scenario driven_cavity_mpiDomainSplitting --solver mpi --numproc 5 4 3 --output false
echo 4 4 3
mpirun -np 48 ./cfd3d --scenario driven_cavity_mpiDomainSplitting --solver mpi --numproc 4 4 3 --output false
