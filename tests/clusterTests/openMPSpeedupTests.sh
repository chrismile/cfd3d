#!/bin/bash

#SBATCH --cluster=mpp2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-type=end
#SBATCH --mail-user=stefan.haas@drei.at
#SBATCH --time=3:00:00
#SBATCH --output=openMPSpeedupTests_%A.log

source /home/hpc/t1221/lu26xag/modules
cd /home/hpc/t1221/lu26xag/cfd3d/build/

#echo OMP_NUM_THREADS=1
#export OMP_NUM_THREADS=1
#./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cpp --linsolver jacobi --output false
#
#echo OMP_NUM_THREADS=2
#export OMP_NUM_THREADS=2
#./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cpp --linsolver jacobi --output false

echo OMP_NUM_THREADS=4
export OMP_NUM_THREADS=4
./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cpp --linsolver jacobi --output false

#echo OMP_NUM_THREADS=8
#export OMP_NUM_THREADS=8
#./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cpp --linsolver jacobi --output false
#
#echo OMP_NUM_THREADS=16
#export OMP_NUM_THREADS=16
#./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cpp --linsolver jacobi --output false
#
#echo OMP_NUM_THREADS=28
#export OMP_NUM_THREADS=28
#./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cpp --linsolver jacobi --output false



