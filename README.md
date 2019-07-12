# cfd3d

A CFD solver using the incompressible Navier-Stokes equations in 3D.
The solver is able to use either OpenMP, MPI or CUDA.

## Building the programm

On Ubuntu 18.04 for example, you can install all necessary packages with this command:

```
sudo apt-get install git cmake libboost-filesystem-dev libnetcdf-dev netcdf-bin libglm-dev libopenmpi-dev
```

For CUDA support, additionally the following packages are necessary.

```
sudo apt-get install nvidia-cuda-dev nvidia-cuda-toolkit
```

Then, to build the program, execute the following commands in the repository directory.

```
mkdir build
cd build
cmake ..
make -j 4
```

This builds the program with only the OpenMP solver. If the user wishes to also build the program with MPI or CUDA
support, run one of the following commands instead.


For MPI support:

```
cmake .. -DUSE_MPI=ON
```

For CUDA support:

```
cmake .. -DUSE_CUDA=ON
```

The program can also be built with both solvers enabled at the same time.

## Running the program

To start the program, the command below can be used.

```
./cfd3d --scenario <scenario-name> --solver <solver-name>
```

The scenario name is the name of one of the scenario files in the 'scenario/' folder (without the file ending).
Here are some examples of how to call the program.

```
./cfd3d --scenario driven_cavity --solver cpp
mpirun -np 8 ./cfd3d --scenario driven_cavity --solver mpi --numproc 2 2 2
./cfd3d --scenario driven_cavity --solver cuda --tracestreamlines true --numparticles 1000
./cfd3d --scenario natural_convection --solver cuda --tracestreaklines true --numparticles 4000
./cfd3d --scenario natural_convection --solver cuda --tracepathlines true
```

The solver name is either 'cpp' for the C++ OpenMP-accelerated solver, 'cuda' for the NVIDIA CUDA solver and 'mpi' for
the MPI solver.
Please note that for the CUDA solver and the MPI solver, the program needs to be built with the necessary flags.

The valid values for all possible arguments are:
* scenario: driven_cavity, flow_over_step, natural_convection, rayleigh_benard_convection_8-2-1,
rayleigh_benard_convection_8-2-2, rayleigh_benard_convection_8-2-4, rayleigh_benard_convection_2d,
single_tower, terrain_1
* solver: cpp, mpi, cuda, sycl
* outputformat: netcdf, vtk (= vtk-binary), vtk-binary, vtk-ascii
* output: true, false (whether to write an output file)
* linsolver: jacobi, sor, gauss-seidel
* tracestreamlines: false, true
* tracestreaklines: false, true
* tracepathlines: false, true
* numparticles: any positive integer number

The standard values for the arguments are:
* scenario: driven_cavity
* solver: cpp
* outputformat: vtk
* outputformat: true
* linsolver: jacobi
* tracestreamlines: false
* tracestreaklines: false
* tracepathlines: false
* numparticles: 1000

Additionally, for the MPI solver, the user MUST also specify the number of processes in x, y and z direction (which must
match the total number of MPI processes):
* numproc: integer x integer x integer

For the CUDA solver, the user CAN also specify the block size in x, y and z direction:
* blocksize: integer x integer x integer

The standard values for the CUDA solver are:
* blocksize: 8 x 8 x 4
