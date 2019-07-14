#!/bin/bash

#SBATCH --cluster=mpp2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-type=end
#SBATCH --mail-user=stefan.haas@drei.at
#SBATCH --time=05:00:00
#SBATCH --output=jacobiSorOpenMPTests_%A.log

source /home/hpc/t1221/lu26xag/modules
cd /home/hpc/t1221/lu26xag/cfd3d/build/
export OMP_NUM_THREADS=28

echo jacobi
./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cpp --linsolver jacobi --output false
echo sor
./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cpp --linsolver sor --output false
echo gauss-seidel
./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cpp --linsolver gauss-seidel --output false
