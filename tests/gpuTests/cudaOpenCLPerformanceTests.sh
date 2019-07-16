#!/bin/bash

# Run: ./cudaOpenCLPerformanceTests.sh > cudaOpenCLPerformanceTests.txt

cd ../../build-float/

echo cuda-rayleigh-benard-float
./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cuda --output false
echo opencl-rayleigh-benard-float
./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver opencl --output false
echo cpp-rayleigh-benard-float
OMP_NUM_THREADS=6 ./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cpp --output false
echo mpi-rayleigh-benard-float
./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver mpi --numproc 3 2 1 --output false

cd ../build-double/

echo cuda-rayleigh-benard-double
./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cuda --output false
echo opencl-rayleigh-benard-double
./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver opencl --output false
echo cpp-rayleigh-benard-double
OMP_NUM_THREADS=6 ./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cpp --output false
echo mpi-rayleigh-benard-double
./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver mpi --numproc 3 2 1 --output false

