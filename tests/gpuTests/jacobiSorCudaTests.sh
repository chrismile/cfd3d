#!/bin/bash

cd /../../build/

echo jacobi
./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cuda --linsolver jacobi --output false
echo sor
./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cuda --linsolver sor --output false
echo gauss-seidel
./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cuda --linsolver gauss-seidel --output false
