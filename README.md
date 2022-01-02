# cfd3d

A CFD solver using the incompressible Navier-Stokes equations in 3D.
The solver is able to use either OpenMP, MPI, CUDA or OpenCL.

Used literature for the theoretical background:
- M. Griebel, T. Dornseifer, and T. Neunhoeffer. SIAM, Philadelphia, 1998. Numerical Simulation in Fluid Dynamics,
a Practical Introduction.


## Building the programm

On Ubuntu 20.04 for example, you can install all necessary packages with this command:

```
sudo apt-get install git cmake libnetcdf-dev libglm-dev libopenmpi-dev
```

For CUDA support, additionally the following packages are necessary.

```
sudo apt-get install nvidia-cuda-dev nvidia-cuda-toolkit
```

For OpenCL support, additionally the following packages are necessary.

```
sudo apt-get install ocl-icd-opencl-dev opencl-headers clinfo
```

The application was tested using NVIDIA's OpenCL implementation, Intel NEO and POCL. However, it should also run with
different implementations. Depending on the OpenCL implementation you want to use, install one of the following packages.


```
sudo apt-get install nvidia-opencl-dev intel-opencl-icd beignet-opencl-icd pocl-opencl-icd
```

intel-opencl-icd is the new OpenCL implementation of Intel that was added in Ubuntu 19.04 to the distribution
repositories. For previous distributions, please use Beignet. POCL is an OpenCL implementation using the user's CPU.


Then, to build the program, execute the following commands in the repository directory.

```
mkdir build
cd build
cmake ..
make -j 4
```

The program will try to automatically find the installed CUDA, OpenCL and MPI library locations.
The program can be built with all solvers enabled at the same time.


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
./cfd3d --scenario driven_cavity --solver cuda
./cfd3d --scenario rayleigh_benard_convection_8-2-1 --solver cuda --tracestreamlines true --numparticles 500
```

The solver name is either 'cpp' for the C++ OpenMP-accelerated solver, 'mpi' for the MPI solver, 'cuda' for the NVIDIA
CUDA solver, and 'opencl' for the OpenCL solver.
Please note that for the MPI solver, CUDA solver and the OpenCL solver, the program needs to be built with the necessary
flags.

The valid values for all possible arguments are:
* scenario: driven_cavity, flow_over_step, natural_convection, rayleigh_benard_convection_8-2-1,
rayleigh_benard_convection_8-2-2, rayleigh_benard_convection_8-2-4,
single_tower, terrain_1, fuji_san, zugspitze
* solver: cpp, mpi, cuda, opencl
* outputformat: netcdf, vtk (= vtk-binary), vtk-binary, vtk-ascii
* output: true, false (whether to write an output file)
* linsolver: jacobi, sor, gauss-seidel
* tracestreamlines: false, true
* numparticles: any positive integer number

The standard values for the arguments are:
* scenario: driven_cavity
* solver: cpp
* outputformat: vtk
* output: true
* linsolver: jacobi
* tracestreamlines: false
* numparticles: 500

Additionally, for the MPI solver, the user MUST also specify the number of processes in x, y and z direction (which must
match the total number of MPI processes):
* numproc: integer integer integer

For the CUDA and OpenCL solvers, the user CAN also specify the block size in x, y and z direction:
* blocksize: integer integer integer

The standard values for the solvers are:
* blocksize: 8 8 4

For the OpenCL solver, additonally to the block size, the user can also specify the ID of the OpenCL platform to use.
If the ID is not specified, it is set to zero. Which platform corresponds to which ID can be found outwith the command
line tool 'clinfo'.
* platformid: integer

For the CUDA and OpenCL solvers, the user CAN also specify the index of the GPU to use for computations (standard is 0):
* gpu: integer
