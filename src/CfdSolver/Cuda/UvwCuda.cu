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

#include <algorithm>
#include <iostream>
#include "UvwCuda.hpp"
#include "CudaDefines.hpp"

__global__ void calculateFghCudaKernel(
        Real Re, Real GX, Real GY, Real GZ, Real alpha, Real beta,
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *T, Real *F, Real *G, Real *H, FlagType *Flag) {
    int i = blockIdx.z * blockDim.z + threadIdx.z + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.x * blockDim.x + threadIdx.x + 1;

    Real d2u_dx2,d2u_dy2,d2u_dz2,
            d2v_dx2,d2v_dy2,d2v_dz2,
            d2w_dx2,d2w_dy2,d2w_dz2;
        
    Real du2_dx,duv_dy,duw_dz,
            duv_dx,dv2_dy,dvw_dz,
            duw_dx,dvw_dy,dw2_dz;
        
    Real Dx = 1/dx, Dy = 1/dy, Dz = 1/dz;
    
    if (i <= imax - 1 && j <= jmax && k <= kmax) {
        if(isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i+1,j,k)])){
            d2u_dx2 = (U[IDXU(i+1,j,k)] - 2*U[IDXU(i,j,k)] + U[IDXU(i-1,j,k)])/(dx*dx);
            d2u_dy2 = (U[IDXU(i,j+1,k)] - 2*U[IDXU(i,j,k)] + U[IDXU(i,j-1,k)])/(dy*dy);
            d2u_dz2 = (U[IDXU(i,j,k+1)] - 2*U[IDXU(i,j,k)] + U[IDXU(i,j,k-1)])/(dz*dz);

            du2_dx = Real(0.25)*Dx*(
                    (U[IDXU(i,j,k)]+U[IDXU(i+1,j,k)])*(U[IDXU(i,j,k)]+U[IDXU(i+1,j,k)]) -
                    (U[IDXU(i-1,j,k)]+U[IDXU(i,j,k)])*(U[IDXU(i-1,j,k)]+U[IDXU(i,j,k)]) +
                    alpha*(
                            (std::abs(U[IDXU(i,j,k)]+U[IDXU(i+1,j,k)])*(U[IDXU(i,j,k)]-U[IDXU(i+1,j,k)]))-
                            (std::abs(U[IDXU(i-1,j,k)]+U[IDXU(i,j,k)])*(U[IDXU(i-1,j,k)]-U[IDXU(i,j,k)]))
                    )
            );

            duv_dy = Real(0.25)*Dy*(
                    (V[IDXV(i,j,k)]+V[IDXV(i+1,j,k)])*(U[IDXU(i,j,k)]+U[IDXU(i,j+1,k)]) -
                    (V[IDXV(i,j-1,k)]+V[IDXV(i+1,j-1,k)])*(U[IDXU(i,j-1,k)]+U[IDXU(i,j,k)]) +
                    alpha*(
                            (std::abs(V[IDXV(i,j,k)]+V[IDXV(i+1,j,k)])*(U[IDXU(i,j,k)]-U[IDXU(i,j+1,k)]))-
                            (std::abs(V[IDXV(i,j-1,k)]+V[IDXV(i+1,j-1,k)])*(U[IDXU(i,j-1,k)]-U[IDXU(i,j,k)]))
                    )
            );

            duw_dz = Real(0.25)*Dz*(
                    (W[IDXW(i,j,k)]+W[IDXW(i+1,j,k)])*(U[IDXU(i,j,k)]+U[IDXU(i,j,k+1)]) -
                    (W[IDXW(i,j,k-1)]+W[IDXW(i+1,j,k-1)])*(U[IDXU(i,j,k-1)]+U[IDXU(i,j,k)]) +
                    alpha*(
                            (std::abs(W[IDXW(i,j,k)]+W[IDXW(i+1,j,k)])*(U[IDXU(i,j,k)]-U[IDXU(i,j,k+1)]))-
                            (std::abs(W[IDXW(i,j,k-1)]+W[IDXW(i+1,j,k-1)])*(U[IDXU(i,j,k-1)]-U[IDXU(i,j,k)]))
                    )
            );

            F[IDXF(i,j,k)] = U[IDXU(i,j,k)] + dt * (
                    (1/Re)*(d2u_dx2+d2u_dy2+d2u_dz2)-
                    du2_dx-duv_dy-duw_dz+
                    GX-(beta/2)*(T[IDXT(i,j,k)]+T[IDXT(i+1,j,k)])*GX
            );
        } else if (B_L(Flag[IDXFLAG(i,j,k)])) {
                F[IDXF(i-1,j,k)] = U[IDXU(i-1,j,k)];
        } else if (B_R(Flag[IDXFLAG(i,j,k)])) {
                F[IDXF(i,j,k)] = U[IDXU(i,j,k)];
        }
    }

    if (i <= imax && j <= jmax - 1 && k <= kmax){
        if(isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i,j+1,k)])){
            d2v_dx2 = (V[IDXV(i+1,j,k)] - 2*V[IDXV(i,j,k)] + V[IDXV(i-1,j,k)])/(dx*dx);
            d2v_dy2 = (V[IDXV(i,j+1,k)] - 2*V[IDXV(i,j,k)] + V[IDXV(i,j-1,k)])/(dy*dy);
            d2v_dz2 = (V[IDXV(i,j,k+1)] - 2*V[IDXV(i,j,k)] + V[IDXV(i,j,k-1)])/(dz*dz);

            duv_dx = Real(0.25)*Dx*(
                    (U[IDXU(i,j,k)]+U[IDXU(i,j+1,k)])*(V[IDXV(i,j,k)]+V[IDXV(i+1,j,k)]) -
                    (U[IDXU(i-1,j,k)]+U[IDXU(i-1,j+1,k)])*(V[IDXV(i-1,j,k)]+V[IDXV(i,j,k)]) +
                    alpha*(
                            (std::abs(U[IDXU(i,j,k)]+U[IDXU(i,j+1,k)])*(V[IDXV(i,j,k)]-V[IDXV(i+1,j,k)]))-
                            (std::abs(U[IDXU(i-1,j,k)]+U[IDXU(i-1,j+1,k)])*(V[IDXV(i-1,j,k)]-V[IDXV(i,j,k)]))
                    )
            );

            dv2_dy = Real(0.25)*Dy*(
                    (V[IDXV(i,j,k)]+V[IDXV(i,j+1,k)])*(V[IDXV(i,j,k)]+V[IDXV(i,j+1,k)]) -
                    (V[IDXV(i,j-1,k)]+V[IDXV(i,j,k)])*(V[IDXV(i,j-1,k)]+V[IDXV(i,j,k)]) +
                    alpha*(
                            (std::abs(V[IDXV(i,j,k)]+V[IDXV(i,j+1,k)])*(V[IDXV(i,j,k)]-V[IDXV(i,j+1,k)]))-
                            (std::abs(V[IDXV(i,j-1,k)]+V[IDXV(i,j,k)])*(V[IDXV(i,j-1,k)]-V[IDXV(i,j,k)]))
                    )
            );

            dvw_dz = Real(0.25)*Dz*(
                    (W[IDXW(i,j,k)]+W[IDXW(i,j+1,k)])*(V[IDXV(i,j,k)]+V[IDXV(i,j,k+1)]) -
                    (W[IDXW(i,j,k-1)]+W[IDXW(i,j+1,k-1)])*(V[IDXV(i,j,k-1)]+V[IDXV(i,j,k)]) +
                    alpha*(
                            (std::abs(W[IDXW(i,j,k)]+W[IDXW(i,j+1,k)])*(V[IDXV(i,j,k)]-V[IDXV(i,j,k+1)]))-
                            (std::abs(W[IDXW(i,j,k-1)]+W[IDXW(i,j+1,k-1)])*(V[IDXV(i,j,k-1)]-V[IDXV(i,j,k)]))
                    )
            );

            G[IDXG(i,j,k)] = V[IDXV(i,j,k)] + dt * (
                    (1/Re)*(d2v_dx2+d2v_dy2+d2v_dz2)-
                    duv_dx-dv2_dy-dvw_dz+
                    GY-(beta/2)*(T[IDXT(i,j,k)]+T[IDXT(i,j+1,k)])*GY
            );
        } else if (B_D(Flag[IDXFLAG(i,j,k)])) {
            G[IDXG(i,j-1,k)] = V[IDXV(i,j-1,k)];
        } else if (B_U(Flag[IDXFLAG(i,j,k)])) {
            G[IDXG(i,j,k)] = V[IDXV(i,j,k)];
        }        
    }

    if (i <= imax && j <= jmax && k <= kmax - 1) {
        if (isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i,j,k+1)])){
            d2w_dx2 = (W[IDXW(i+1,j,k)] - 2*W[IDXW(i,j,k)] + W[IDXW(i-1,j,k)])/(dx*dx);
            d2w_dy2 = (W[IDXW(i,j+1,k)] - 2*W[IDXW(i,j,k)] + W[IDXW(i,j-1,k)])/(dy*dy);
            d2w_dz2 = (W[IDXW(i,j,k+1)] - 2*W[IDXW(i,j,k)] + W[IDXW(i,j,k-1)])/(dz*dz);

            duw_dx = Real(0.25)*Dx*(
                    (U[IDXU(i,j,k)]+U[IDXU(i,j,k+1)])*(W[IDXW(i,j,k)]+W[IDXW(i+1,j,k)]) -
                    (U[IDXU(i-1,j,k)]+U[IDXU(i-1,j,k+1)])*(W[IDXW(i-1,j,k)]+W[IDXW(i,j,k)]) +
                    alpha*(
                            (std::abs(U[IDXU(i,j,k)]+U[IDXU(i,j,k+1)])*(W[IDXW(i,j,k)]-W[IDXW(i+1,j,k)]))-
                            (std::abs(U[IDXU(i-1,j,k)]+U[IDXU(i-1,j,k+1)])*(W[IDXW(i-1,j,k)]-W[IDXW(i,j,k)]))
                    )
            );

            dvw_dy = Real(0.25)*Dy*(
                    (V[IDXV(i,j,k)]+V[IDXV(i,j,k+1)])*(W[IDXW(i,j,k)]+W[IDXW(i,j+1,k)]) -
                    (V[IDXV(i,j-1,k)]+V[IDXV(i,j-1,k+1)])*(W[IDXW(i,j-1,k)]+W[IDXW(i,j,k)]) +
                    alpha*(
                            (std::abs(V[IDXV(i,j,k)]+V[IDXV(i,j,k+1)])*(W[IDXW(i,j,k)]-W[IDXW(i,j+1,k)]))-
                            (std::abs(V[IDXV(i,j-1,k)]+V[IDXV(i,j-1,k+1)])*(W[IDXW(i,j-1,k)]-W[IDXW(i,j,k)]))
                    )
            );

            dw2_dz = Real(0.25)*Dz*(
                    (W[IDXW(i,j,k)]+W[IDXW(i,j,k+1)])*(W[IDXW(i,j,k)]+W[IDXW(i,j,k+1)]) -
                    (W[IDXW(i,j,k-1)]+W[IDXW(i,j,k)])*(W[IDXW(i,j,k-1)]+W[IDXW(i,j,k)]) +
                    alpha*(
                            (std::abs(W[IDXW(i,j,k)]+W[IDXW(i,j,k+1)])*(W[IDXW(i,j,k)]-W[IDXW(i,j,k+1)]))-
                            (std::abs(W[IDXW(i,j,k-1)]+W[IDXW(i,j,k)])*(W[IDXW(i,j,k-1)]-W[IDXW(i,j,k)]))
                    )
            );

            H[IDXH(i,j,k)] = W[IDXW(i,j,k)] +  dt * (
                    (1/Re)*(d2w_dx2+d2w_dy2+d2w_dz2)-
                    duw_dx-dvw_dy-dw2_dz+
                    GZ-(beta/2)*(T[IDXT(i,j,k)]+T[IDXT(i,j,k+1)])*GZ
            );
        } else if (B_B(Flag[IDXFLAG(i,j,k)])) {
            H[IDXH(i,j,k-1)] = W[IDXW(i,j,k-1)];
        } else if (B_F(Flag[IDXFLAG(i,j,k)])) {
            H[IDXH(i,j,k)] = W[IDXW(i,j,k)];
        }    
    }
}

__global__ void setFBoundariesCuda(int imax, int jmax, int kmax, Real *U, Real *F) {

    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.x * blockDim.x + threadIdx.x + 1;

    // Set the boundary values for F on the y-z-planes.
    if (j <= jmax && k <= kmax){
        F[IDXF(0,j,k)] = U[IDXU(0,j,k)];         
        F[IDXF(imax,j,k)] = U[IDXU(imax,j,k)];
    }
}

__global__ void setGBoundariesCuda(int imax, int jmax, int kmax, Real *V, Real *G) {

    int i = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.x * blockDim.x + threadIdx.x + 1;

    // Set the boundary values for G on the x-z-planes.
    if (i <= imax && k <= kmax){
        G[IDXG(i,0,k)] = V[IDXV(i,0,k)];         
        G[IDXG(i,jmax,k)] = V[IDXV(i,jmax,k)];
    }
}

__global__ void setHBoundariesCuda(int imax, int jmax, int kmax, Real *W, Real *H) {

    int i = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int j = blockIdx.x * blockDim.x + threadIdx.x + 1;

    // Set the boundary values for G on the x-z-planes.
    if (i <= imax && j <= jmax){
        H[IDXH(i,j,0)] = W[IDXW(i,j,0)];         
        H[IDXH(i,j,kmax)] = W[IDXW(i,j,kmax)];
    }
}


void calculateFghCuda(
        Real Re, Real GX, Real GY, Real GZ, Real alpha, Real beta,
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        int blockSizeX, int blockSizeY, int blockSizeZ,
        Real *U, Real *V, Real *W, Real *T, Real *F, Real *G, Real *H, FlagType *Flag) {

    dim3 dimBlock(blockSizeX, blockSizeY, blockSizeZ);
    dim3 dimGrid(iceil(kmax,dimBlock.x),iceil(jmax,dimBlock.y),iceil(imax,dimBlock.z));
    calculateFghCudaKernel<<<dimGrid,dimBlock>>>(
            Re, GX, GY, GZ, alpha, beta, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, T, F, G, H, Flag);

    dim3 dimBlock2D(blockSizeX, blockSizeY);
    dim3 dimGrid_y_z(iceil(kmax,dimBlock2D.x),iceil(jmax,dimBlock2D.y));
    setFBoundariesCuda<<<dimGrid_y_z,dimBlock2D>>>(imax, jmax, kmax, U, F);

    dim3 dimGrid_x_z(iceil(kmax,dimBlock2D.x),iceil(imax,dimBlock2D.y));
    setGBoundariesCuda<<<dimGrid_x_z,dimBlock2D>>>(imax, jmax, kmax, V, G);

    dim3 dimGrid_x_y(iceil(jmax,dimBlock2D.x),iceil(imax,dimBlock2D.y));
    setHBoundariesCuda<<<dimGrid_x_y,dimBlock2D>>>(imax, jmax, kmax, W, H);
        
}

__global__ void calculateRsCuda(
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *F, Real *G, Real *H, Real *RS) {
    int i = blockIdx.z * blockDim.z + threadIdx.z + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.x * blockDim.x + threadIdx.x + 1;

    if (i <= imax && j <= jmax && k <= kmax){
        RS[IDXRS(i, j, k)] = ((F[IDXF(i, j, k)] - F[IDXF(i - 1, j, k)]) / dx +
                (G[IDXG(i, j, k)] - G[IDXG(i, j - 1, k)]) / dy +
                (H[IDXH(i, j, k)] - H[IDXH(i, j, k - 1)]) / dz) / dt;
    }
}

/**
 * Reference: Based on kernel 4 from https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
 * @param input The array of input values (of size 'sizeOfInput').
 * @param output The output array (of size iceil(numberOfBlocksI, blockSize1D*2)).
 * @param sizeOfInput The number of input values.
 */
__global__ void calculateMaximum(Real *input, Real *output, int sizeOfInput, int blockSize1D) {
    extern __shared__ Real sdata[];

    unsigned int threadID = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x * 2 + threadIdx.x;

    // Copy the data to the shared memory and do the first reduction step.
    if (i + blockDim.x < sizeOfInput){
#ifdef REAL_DOUBLE
    sdata[threadID] = fmax(fabs(input[i]), fabs(input[i + blockDim.x]));
#else
    sdata[threadID] = fmaxf(fabsf(input[i]), fabsf(input[i + blockDim.x]));
#endif
    } else if (i < sizeOfInput){
        sdata[threadID] = input[i];
    } else{
        sdata[threadID] = 0;
    }
    __syncthreads();

    // Do the reduction in the shared memory.
    for (unsigned int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadID < stride) {
#ifdef REAL_DOUBLE
            sdata[threadID] = fmax(fabs(sdata[threadID]), fabs(sdata[threadID + stride]));
#else
            sdata[threadID] = fmaxf(fabsf(sdata[threadID]), fabsf(sdata[threadID + stride]));
#endif
        }
        __syncthreads();
    }

    // Write the result for this block to global memory.
    if (threadID == 0) {
        output[blockIdx.x] = sdata[0];
    }
}

void calculateDtCuda(
        Real Re, Real Pr, Real tau,
        Real &dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax, int blockSize1D,
        Real *U, Real *V, Real *W,
        Real *cudaReductionArrayU1, Real *cudaReductionArrayU2,
        Real *cudaReductionArrayV1, Real *cudaReductionArrayV2,
        Real *cudaReductionArrayW1, Real *cudaReductionArrayW2,
        bool useTemperature) {
    Real uMaxAbs = Real(0.0), vMaxAbs = Real(0.0), wMaxAbs = Real(0.0);
    
    Real *U_reductionInput = U;
    Real *V_reductionInput = V;
    Real *W_reductionInput = W;
    Real *U_reductionOutput = cudaReductionArrayU1;
    Real *V_reductionOutput = cudaReductionArrayV1;
    Real *W_reductionOutput = cudaReductionArrayW1;

    int numberOfBlocksI = (imax+1)*(jmax+2)*(kmax+2);
    int numberOfBlocksJ = (imax+2)*(jmax+1)*(kmax+2);
    int numberOfBlocksK = (imax+2)*(jmax+2)*(kmax+1);
    int inputSizeI;
    int inputSizeJ;
    int inputSizeK;
    bool finished = false;

    //const int blockSize1D = 256;

    int iteration = 0;
    while (!finished) {
        inputSizeI = numberOfBlocksI;
        inputSizeJ = numberOfBlocksJ;
        inputSizeK = numberOfBlocksK;
        numberOfBlocksI = iceil(numberOfBlocksI, blockSize1D*2);
        numberOfBlocksJ = iceil(numberOfBlocksJ, blockSize1D*2);
        numberOfBlocksK = iceil(numberOfBlocksK, blockSize1D*2);

        if (inputSizeI != 1) {
            int sharedMemorySize = blockSize1D * sizeof(Real);
            calculateMaximum<<<numberOfBlocksI, blockSize1D, sharedMemorySize>>>(
                    U_reductionInput, U_reductionOutput, inputSizeI, blockSize1D);
            if (iteration % 2 == 0) {
                U_reductionInput = cudaReductionArrayU1;
                U_reductionOutput = cudaReductionArrayU2;
            } else {
                U_reductionInput = cudaReductionArrayU2;
                U_reductionOutput = cudaReductionArrayU1;
            }
        }
        if (inputSizeJ != 1) {
            int sharedMemorySize = blockSize1D * sizeof(Real);
            calculateMaximum <<<numberOfBlocksJ, blockSize1D, sharedMemorySize>>> (
                    V_reductionInput, V_reductionOutput, inputSizeJ, blockSize1D);
            if (iteration % 2 == 0) {
                V_reductionInput = cudaReductionArrayV1;
                V_reductionOutput = cudaReductionArrayV2;
            } else {
                V_reductionInput = cudaReductionArrayV2;
                V_reductionOutput = cudaReductionArrayV1;
            }
        }
        if (inputSizeK != 1) {
            int sharedMemorySize = blockSize1D * sizeof(Real);
            calculateMaximum <<<numberOfBlocksK, blockSize1D, sharedMemorySize>>> (
                    W_reductionInput, W_reductionOutput, inputSizeK, blockSize1D);
            if (iteration % 2 == 0) {
                W_reductionInput = cudaReductionArrayW1;
                W_reductionOutput = cudaReductionArrayW2;
            } else {
                W_reductionInput = cudaReductionArrayW2;
                W_reductionOutput = cudaReductionArrayW1;
            }
        }


        if (numberOfBlocksI == 1 && numberOfBlocksJ == 1 && numberOfBlocksK == 1) {
            finished = true;
        }
        iteration++;
    }

    cudaMemcpy(&uMaxAbs, U_reductionInput, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(&vMaxAbs, V_reductionInput, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(&wMaxAbs, W_reductionInput, sizeof(Real), cudaMemcpyDeviceToHost);

    if (tau < Real(0.0)) {
        // Constant time step manually specified in configuration file. Check for stability.
        assert(2 / Re * dt < dx * dx * dy * dy * dz * dz / (dx * dx + dy * dy + dz * dz));
        assert(uMaxAbs * dt < dx);
        assert(vMaxAbs * dt < dy);
        assert(wMaxAbs * dt < dz);
        if (useTemperature){
            assert(dt < (Re*Pr/2)*(1/((1/(dx*dx))+1/(dy*dy)+1/(dz*dz))));
        }
        return;
    }

    dt = std::min(dx / uMaxAbs, dy / vMaxAbs);
    dt = std::min(dt, dz / wMaxAbs);
    dt = std::min(dt, (Re / Real(2.0)) * (Real(1.0) / (Real(1.0) / (dx*dx) + Real(1.0) / (dy*dy)
            + Real(1.0) / (dz*dz))));
    if (useTemperature){
        dt = std::min(dt, (Re * Pr / Real(2.0)) * (Real(1.0) / (Real(1.0) / (dx*dx) + Real(1.0) / (dy*dy)
                + Real(1.0) / (dz*dz))));
    }
    dt = tau * dt;
}

__global__ void calculateUvwCuda(
    Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
    Real *U, Real *V, Real *W, Real *F, Real *G, Real *H, Real *P, FlagType *Flag) {

    int i = blockIdx.z * blockDim.z + threadIdx.z + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.x * blockDim.x + threadIdx.x + 1;

    if (i <= imax - 1 && j <= jmax && k <= kmax){
        if(isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i+1,j,k)])){
            U[IDXU(i, j, k)] = F[IDXF(i, j, k)] - dt / dx * (P[IDXP(i + 1, j, k)] - P[IDXP(i, j, k)]);
        }
    }

    if (i <= imax && j <= jmax - 1 && k <= kmax){
        if(isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i,j+1,k)])){
            V[IDXV(i, j, k)] = G[IDXG(i, j, k)] - dt / dy * (P[IDXP(i, j + 1, k)] - P[IDXP(i, j, k)]);
        }
    }

    if (i <= imax && j <= jmax && k <= kmax - 1){
        if(isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i,j,k+1)])){
            W[IDXW(i, j, k)] = H[IDXH(i, j, k)] - dt / dz * (P[IDXP(i, j, k + 1)] - P[IDXP(i, j, k)]);
        }
    }
}

__global__ void calculateTemperatureCuda(
    Real Re, Real Pr, Real alpha,
    Real dt, Real dx, Real dy, Real dz,
    int imax, int jmax, int kmax,
    Real *U, Real *V, Real *W, Real *T, Real *T_temp, FlagType *Flag) {

    int i = blockIdx.z * blockDim.z + threadIdx.z + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.x * blockDim.x + threadIdx.x + 1;

    Real duT_dx, dvT_dy, dwT_dz, d2T_dx2, d2T_dy2, d2T_dz2;

    if (i <= imax && j <= jmax && k <= kmax){
        if(isFluid(Flag[IDXFLAG(i,j,k)])){
            duT_dx = 1 / dx * (
                    U[IDXU(i, j, k)] * ((T_temp[IDXT(i, j, k)] + T_temp[IDXT(i + 1, j, k)]) / 2) -
                    U[IDXU(i - 1, j, k)] * ((T_temp[IDXT(i - 1, j, k)] + T_temp[IDXT(i, j, k)]) / 2) +
                    alpha * (
                            std::abs(U[IDXU(i, j, k)])*((T_temp[IDXT(i, j, k)] - T_temp[IDXT(i + 1, j, k)]) / 2) -
                            std::abs(U[IDXU(i - 1, j, k)])*((T_temp[IDXT(i - 1, j, k)] - T_temp[IDXT(i, j, k)]) / 2)
                    )
            );

            dvT_dy = 1 / dy * (
                    V[IDXV(i, j, k)] * ((T_temp[IDXT(i, j, k)] + T_temp[IDXT(i, j + 1, k)]) / 2) -
                    V[IDXV(i, j - 1, k)] * ((T_temp[IDXT(i, j - 1, k)] + T_temp[IDXT(i, j, k)]) / 2) +
                    alpha * (
                            std::abs(V[IDXV(i, j, k)])*((T_temp[IDXT(i, j, k)] - T_temp[IDXT(i, j + 1, k)]) / 2) -
                            std::abs(V[IDXV(i, j - 1, k)])*((T_temp[IDXT(i, j - 1, k)] - T_temp[IDXT(i, j, k)]) / 2)
                    )
            );

            dwT_dz = 1 / dz * (
                    W[IDXW(i, j, k)] * ((T_temp[IDXT(i, j, k)] + T_temp[IDXT(i, j, k + 1)]) / 2) -
                    W[IDXW(i, j, k - 1)] * ((T_temp[IDXT(i, j, k - 1)] + T_temp[IDXT(i, j, k)]) / 2) +
                    alpha * (
                            std::abs(W[IDXW(i, j, k)])*((T_temp[IDXT(i, j, k)] - T_temp[IDXT(i, j, k + 1)]) / 2) -
                            std::abs(W[IDXW(i, j, k - 1)])*((T_temp[IDXT(i, j, k - 1)] - T_temp[IDXT(i, j, k)]) / 2)
                    )
            );

            d2T_dx2 =
                    (T_temp[IDXT(i + 1, j, k)] - 2 * T_temp[IDXT(i, j, k)] + T_temp[IDXT(i - 1, j, k)]) / (dx*dx);

            d2T_dy2 =
                    (T_temp[IDXT(i, j + 1, k)] - 2 * T_temp[IDXT(i, j, k)] + T_temp[IDXT(i, j - 1, k)]) / (dy*dy);

            d2T_dz2 =
                    (T_temp[IDXT(i, j, k + 1)] - 2 * T_temp[IDXT(i, j, k)] + T_temp[IDXT(i, j, k - 1)]) / (dz*dz);

            T[IDXT(i, j, k)] = T_temp[IDXT(i, j, k)] + dt * (
                    (1 / (Re*Pr))*(d2T_dx2 + d2T_dy2 + d2T_dz2) -
                    duT_dx -
                    dvT_dy -
                    dwT_dz
            );
        }
    }
}
