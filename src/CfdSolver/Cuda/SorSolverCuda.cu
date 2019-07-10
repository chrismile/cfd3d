/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2019, Christoph Neuhauser, Stefan Haas, Paul Ng
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <iostream>
#include "SorSolverCuda.hpp"
#include "CudaDefines.hpp"

__global__ void set_x_y_planes_pressure_boundaries(
        int imax, int jmax, int kmax, Real *P) {
    int i = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int j = blockIdx.x * blockDim.x + threadIdx.x + 1;

    // Set the boundary values for the pressure on the x-y-planes.
    if (i <= imax && j <= jmax) {
        P[IDXP(i, j, 0)] = P[IDXP(i, j, 1)];
        P[IDXP(i, j, kmax + 1)] = P[IDXP(i, j, kmax)];
    }
}

__global__ void set_x_z_planes_pressure_boundaries(
        int imax, int jmax, int kmax, Real *P) {
    int i = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.x * blockDim.x + threadIdx.x + 1;
    // Set the boundary values for the pressure on the x-z-planes.
    if (i <= imax && k <= jmax) {
        P[IDXP(i, 0, k)] = P[IDXP(i, 1, k)];
        P[IDXP(i, jmax + 1, k)] = P[IDXP(i, jmax, k)];
    }
}

__global__ void set_y_z_planes_pressure_boundaries(
        int imax, int jmax, int kmax, Real *P) {
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.x * blockDim.x + threadIdx.x + 1;
    // Set the boundary values for the pressure on the y-z-planes.
    if (j <= imax && k <= jmax) {
        P[IDXP(0, j, k)] = P[IDXP(1, j, k)];
        P[IDXP(imax + 1, j, k)] = P[IDXP(imax, j, k)];
    }
}

/**
 * Reference: Based on kernel 4 from https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
 * @param input The array of input values (of size 'sizeOfInput').
 * @param output The output array (of size iceil(numberOfBlocksI, blockSize*blockSize*2)).
 * @param sizeOfInput The number of input values.
 */
template<class T>
__global__ void reduceSumCudaKernel(T *input, T *output, int sizeOfInput) {
    extern __shared__ unsigned int sharedMemory[];
    T* sdata = (T*)sharedMemory;

    unsigned int threadID = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x * 2 + threadIdx.x;

    // Copy the data to the shared memory and do the first reduction step.
    if (i + blockDim.x < sizeOfInput){
        sdata[threadID] = input[i] + input[i + blockDim.x];
    } else if (i < sizeOfInput){
        sdata[threadID] = input[i];
    } else{
        sdata[threadID] = 0;
    }
    __syncthreads();

    // Do the reduction in the shared memory.
    for (unsigned int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadID < stride) {
            sdata[threadID] += sdata[threadID + stride];
        }
        __syncthreads();
    }

    // Write the result for this block to global memory.
    if (threadID == 0) {
        output[blockIdx.x] = sdata[0];
    }
}

/*
 * Overwrites the contents of the passed reduction arrays.
 */
template<class T>
T reduceSumCuda(T *input, unsigned int numValues, T *cudaReductionHelperArray) {
    /*cudaDeviceSynchronize();
    T *arrayCpu = new T[numValues];
    memset(arrayCpu, 0, sizeof(T)*numValues);
    int n = sizeof(T)*numValues;
    cudaMemcpy(arrayCpu, input, n, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();*/


    T sumValue = T(0);

    T *reductionInput = input;
    T *reductionOutput = cudaReductionHelperArray;

    int numberOfBlocks = numValues;
    int inputSize;
    bool finished = false;

    //const int blockSize1D = 256;

    int iteration = 0;
    while (!finished) {
        inputSize = numberOfBlocks;
        numberOfBlocks = iceil(numberOfBlocks, blockSize*blockSize*2);

        if (inputSize != 1) {
            int sharedMemorySize = blockSize*blockSize * sizeof(T);
            reduceSumCudaKernel<<<numberOfBlocks, blockSize*blockSize, sharedMemorySize>>>(
                    reductionInput, reductionOutput, inputSize);
            if (iteration % 2 == 0) {
                reductionInput = cudaReductionHelperArray;
                reductionOutput = input;
            } else {
                reductionInput = input;
                reductionOutput = cudaReductionHelperArray;
            }
        }

        if (numberOfBlocks == 1) {
            finished = true;
        }
        iteration++;
    }

    cudaMemcpy(&sumValue, reductionInput, sizeof(T), cudaMemcpyDeviceToHost);


    /*T sumValueCpu = T(0);

    #pragma omp parallel for reduction(+: sumValueCpu)
    for (int i = 0; i < numValues; i++) {
        sumValueCpu += arrayCpu[i];
    }

    std::cout << "Sum GPU: " << sumValue << std::endl;
    std::cout << "Sum CPU: " << sumValueCpu << std::endl;

    delete[] arrayCpu;*/

    return sumValue;
}


__global__ void sorSolverIterationCuda(
        Real omg, Real dx, Real dy, Real dz, Real coeff, int imax, int jmax, int kmax,
        Real *P, Real *P_temp, Real *RS, FlagType *Flag) {
    int i = blockIdx.z * blockDim.z + threadIdx.z + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.x * blockDim.x + threadIdx.x + 1;

    if (i <= imax && j <= jmax && k <= kmax) {
        P[IDXP(i, j, k)] = (Real(1.0) - omg) * P_temp[IDXP(i, j, k)] + coeff *
                ((P_temp[IDXP(i + 1, j, k)] + P_temp[IDXP(i - 1, j, k)]) / (dx * dx)
                 + (P_temp[IDXP(i, j + 1, k)] + P_temp[IDXP(i, j - 1, k)]) / (dy * dy)
                 + (P_temp[IDXP(i, j, k + 1)] + P_temp[IDXP(i, j, k - 1)]) / (dz * dz)
                 - RS[IDXRS(i, j, k)]);
    }
}

__global__ void sorSolverComputeResidualArrayCuda(
        Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *P, Real *RS, FlagType *Flag, Real *residualArray, unsigned int *numFluidCellsArray) {
    int i = blockIdx.z * blockDim.z + threadIdx.z + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.x * blockDim.x + threadIdx.x + 1;

    if (i <= imax && j <= jmax && k <= kmax){
        int arrayIndex1D = (i-1)*jmax*kmax + (j-1)*kmax + k-1;
        if (isFluid(Flag[IDXFLAG(i,j,k)])){
            residualArray[arrayIndex1D] = SQR(
                    (P[IDXP(i+1,j,k)] - Real(2.0)*P[IDXP(i,j,k)] + P[IDXP(i-1,j,k)])/(dx*dx)
                    + (P[IDXP(i,j+1,k)] - Real(2.0)*P[IDXP(i,j,k)] + P[IDXP(i,j-1,k)])/(dy*dy)
                    + (P[IDXP(i,j,k+1)] - Real(2.0)*P[IDXP(i,j,k)] + P[IDXP(i,j,k-1)])/(dz*dz)
                    - RS[IDXRS(i,j,k)]
            );
            numFluidCellsArray[arrayIndex1D] = 1;
        } else {
            residualArray[arrayIndex1D] = Real(0.0);
            numFluidCellsArray[arrayIndex1D] = 0;
        }
    }
}


void sorSolverCuda(
        Real omg, Real eps, int itermax,
        Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *P, Real *P_temp, Real *RS, FlagType *Flag,
        Real *cudaReductionArrayResidual1, Real *cudaReductionArrayResidual2,
        unsigned int *cudaReductionArrayNumCells1, unsigned int *cudaReductionArrayNumCells2) {
    omg = 1.0;
    const Real coeff = omg / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));

    Real residual = 1e9;
    int it = 0;
    while (it < itermax && residual > eps) {
        dim3 dimBlock2D(blockSize, blockSize);
        dim3 dimGrid_x_y(iceil(jmax, dimBlock2D.x), iceil(imax, dimBlock2D.y));
        set_x_y_planes_pressure_boundaries<<<dimGrid_x_y, dimBlock2D>>>(imax, jmax, kmax, P);

        dim3 dimGrid_x_z(iceil(kmax, dimBlock2D.x), iceil(imax, dimBlock2D.y));
        set_x_z_planes_pressure_boundaries<<<dimGrid_x_z, dimBlock2D>>>(imax, jmax, kmax, P);

        dim3 dimGrid_y_z(iceil(kmax, dimBlock2D.x), iceil(jmax, dimBlock2D.y));
        set_y_z_planes_pressure_boundaries<<<dimGrid_y_z, dimBlock2D>>>(imax, jmax, kmax, P);

        cudaMemcpy(P_temp, P, sizeof(Real) * (imax + 2) * (jmax + 2) * (kmax + 2), cudaMemcpyDeviceToDevice);

        dim3 dimBlock(blockSize, blockSize);
        dim3 dimGrid(iceil(kmax, dimBlock.x), iceil(jmax, dimBlock.y), iceil(imax, dimBlock.z));
        sorSolverIterationCuda<<<dimGrid, dimBlock>>>(
                omg, dx, dy, dz, coeff, imax, jmax, kmax, P, P_temp, RS, Flag);

        sorSolverComputeResidualArrayCuda<<<dimGrid, dimBlock>>>(
                dx, dy, dz, imax, jmax, kmax, P, RS, Flag, cudaReductionArrayResidual1, cudaReductionArrayNumCells1);

        residual = reduceSumCuda(
                cudaReductionArrayResidual1, imax*jmax*kmax, cudaReductionArrayResidual2);
        unsigned int numFluidCells = reduceSumCuda(
                cudaReductionArrayNumCells1, imax*jmax*kmax, cudaReductionArrayNumCells2);
        residual = std::sqrt(residual / Real(numFluidCells));

        it++;
    }

    if ((residual > eps && it == itermax) || std::isnan(residual)) {
        std::cout << "\nSOR solver reached maximum number of iterations without converging (res: "
                  << residual << ")." << std::endl;
    }
    if (std::isnan(residual)) {
        std::cout << "\nResidual in SOR solver is not a number." << std::endl;
        exit(1);
    }
}
