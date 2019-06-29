# cfd3d

A CFD solver using the incompressible Navier-Stokes equations in 3D.
The solver is able to use either OpenMP or CUDA.

## Building the programm

On Ubuntu 18.04 for example, you can install all necessary packages with this command:

```
sudo apt-get install git cmake libboost-filesystem-dev libnetcdf-dev netcdf-bin libglm-dev
```

Then, to build the program, execute the following commands in the repository directory.

```
mkdir build
cd build
cmake ..
make -j 4
```

This builds the program with only the OpenMP solver. If the user wishes to also build the program with CUDA or SYCL
support, run one of the following commands instead.


For CUDA support:

```
cmake .. -DUSE_CUDA=ON
```

For SYCL support using the patched LLVM compiler from Intel (for more details see https://github.com/intel/llvm/blob/sycl/sycl/doc/GetStartedWithSYCLCompiler.md):

```
CXX=$SYCL_HOME/build/bin/clang++ CXXFLAGS="-fsycl" LDFLAGS="-lOpenCL" cmake .. -DUSE_SYCL=ON
```

For SYCL support using hipSYCL with the CUDA backend (for more details see https://github.com/illuhad/hipSYCL):

```
CXX=syclcc-clang CXXFLAGS="--hipsycl-platform=cuda" cmake .. -DUSE_SYCL=ON
```

## Running the program

To start the program, the command below can be used.

```
./cfd3d --scenario <scenario-name> --solver <solver-name> --tracestreamlines <true-or-false> \
--tracestreaklines <true-or-false> --tracepathlines <true-or-false> --numparticles <numparticles>
```

The scenario name is the name of one of the scenario files in the 'scenario/' folder (without the file ending).
Here are some examples how to call the program.

```
./cfd3d --scenario driven_cavity --solver cpp
./cfd3d --scenario driven_cavity --solver cuda --tracestreamlines true --numparticles 1000
./cfd3d --scenario natural_convection --solver cuda --tracestreaklines true --numparticles 4000
./cfd3d --scenario natural_convection --solver cuda --tracepathlines true
```

The solver name is either 'cpp' for the C++ OpenMP-accelerated solver, 'cuda' for the NVIDIA CUDA solver, or 'sycl' for
the SYCL solver. Please note that for the CUDA solver and the SYCL solver, the program needs to be built with the
necessary flags.

The standard values for the arguments are:
* scenario: driven_cavity
* solver: cpp
* tracestreamlines: false
* tracestreaklines: false
* tracepathlines: false
* numparticles: 1000
