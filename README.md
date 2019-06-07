# cfd3d

A CFD solver using the incompressible Navier-Stokes equations in 3D.
The solver is able to use either OpenMP or CUDA.

## Building and running the programm

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
CXX=$SYCL_HOME/build/bin/clang++ CXXFLAGS="-fsycl -lOpenCL" cmake .. -DUSE_SYCL=ON
```

For SYCL support using hipSYCL with the CUDA backend (for more details see https://github.com/illuhad/hipSYCL):

```
CXX=syclcc-clang CXXFLAGS="--hipsycl-platform=cuda" cmake .. -DUSE_SYCL=ON
```
