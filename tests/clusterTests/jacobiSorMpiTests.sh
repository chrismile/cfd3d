#!/bin/bash

#SBATCH --cluster=mpp2
#SBATCH --nodes=1
#SBATCH --ntasks=27
#SBATCH --mail-type=end
#SBATCH --mail-user=stefan.haas@drei.at
#SBATCH --time=05:00:00
#SBATCH --output=jacobiSorMpiTests_%A.log

source /home/hpc/t1221/lu26xag/modules
cd /home/hpc/t1221/lu26xag/cfd3d/build/

echo jacobi
mpirun -np 27 ./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver mpi --numproc 3 3 3 --linsolver jacobi --output false
echo sor
mpirun -np 27 ./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver mpi --numproc 3 3 3 --linsolver sor --output false
echo gauss-seidel
mpirun -np 27 ./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver mpi --numproc 3 3 3 --linsolver gauss-seidel --output false
