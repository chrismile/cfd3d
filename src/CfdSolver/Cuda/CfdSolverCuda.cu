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

#include <cstring>
#include <iostream>
#include "BoundaryValuesCuda.hpp"
#include "UvwCuda.hpp"
#include "SorSolverCuda.hpp"
#include "CfdSolverCuda.hpp"
#include "../../Defines.hpp"

CfdSolverCuda::CfdSolverCuda(int gpuId, int blockSizeX, int blockSizeY, int blockSizeZ, int blockSize1D) {
    this->gpuId = gpuId;
    this->blockSizeX = blockSizeX;
    this->blockSizeY = blockSizeY;
    this->blockSizeZ = blockSizeZ;
    this->blockSize1D = blockSize1D;
}

void CfdSolverCuda::initialize(
        const std::string &scenarioName, LinearSystemSolverType linearSystemSolverType, bool shallWriteOutput,
        Real Re, Real Pr, Real omg, Real eps, int itermax, Real alpha, Real beta, Real dt, Real tau,
        Real GX, Real GY, Real GZ, bool useTemperature, Real T_h, Real T_c,
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz,
        Real *U, Real *V, Real *W, Real *P, Real *T, uint32_t *Flag) {
    this->scenarioName = scenarioName;
    this->linearSystemSolverType = linearSystemSolverType;
    this->shallWriteOutput = shallWriteOutput;
    this->Re = Re;
    this->Pr = Pr;
    this->omg = omg;
    this->eps = eps;
    this->itermax = itermax;
    this->alpha = alpha;
    this->beta = beta;
    this->dt = dt;
    this->tau = tau;
    this->GX = GX;
    this->GY = GY;
    this->GZ = GZ;
    this->useTemperature = useTemperature;
    this->T_h = T_h;
    this->T_c = T_c;
    this->imax = imax;
    this->jmax = jmax;
    this->kmax = kmax;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;

    int numDevices = 0;
    cudaGetDeviceCount(&numDevices);
    if (numDevices == 0) {
        std::cerr << "Fatal error in CfdSolverCuda::initialize: No CUDA devices were found." << std::endl;
        exit(1);
    }
    if (gpuId >= numDevices) {
        std::cerr << "Error in CfdSolverCuda::initialize: Invalid device ID specified. Setting device ID to 0."
                  << std::endl;
        gpuId = 0;
    }
    cudaSetDevice(gpuId);

    // Create all arrays for the simulation.
    cudaMalloc(&this->U, (imax+1)*(jmax+2)*(kmax+2)*sizeof(Real));
    cudaMalloc(&this->V, (imax+2)*(jmax+1)*(kmax+2)*sizeof(Real));
    cudaMalloc(&this->W, (imax+2)*(jmax+2)*(kmax+1)*sizeof(Real));
    cudaMalloc(&this->P, (imax+2)*(jmax+2)*(kmax+2)*sizeof(Real));
    cudaMalloc(&this->P_temp, (imax+2)*(jmax+2)*(kmax+2)*sizeof(Real));
    cudaMalloc(&this->T, (imax+2)*(jmax+2)*(kmax+2)*sizeof(Real));
    cudaMalloc(&this->T_temp, (imax+2)*(jmax+2)*(kmax+2)*sizeof(Real));
    cudaMalloc(&this->F, (imax+1)*(jmax+1)*(kmax+1)*sizeof(Real));
    cudaMalloc(&this->G, (imax+1)*(jmax+1)*(kmax+1)*sizeof(Real));
    cudaMalloc(&this->H, (imax+1)*(jmax+1)*(kmax+1)*sizeof(Real));
    cudaMalloc(&this->RS, (imax+1)*(jmax+1)*(kmax+1)*sizeof(Real));
    cudaMalloc(&this->Flag, (imax+2)*(jmax+2)*(kmax+2)*sizeof(unsigned int));

    int cudaReductionArrayUSize = iceil((imax+1)*(jmax+2)*(kmax+2), blockSize1D*2);
    int cudaReductionArrayVSize = iceil((imax+1)*(jmax+2)*(kmax+2), blockSize1D*2);
    int cudaReductionArrayWSize = iceil((imax+1)*(jmax+2)*(kmax+2), blockSize1D*2);
    int cudaReductionArrayResidualSize1 = iceil(imax*jmax*kmax, blockSize1D*2)*blockSize1D*2;
    int cudaReductionArrayResidualSize2 = iceil(imax*jmax*kmax, blockSize1D*2);
    cudaMalloc(&cudaReductionArrayU1, cudaReductionArrayUSize*sizeof(Real));
    cudaMalloc(&cudaReductionArrayU2, cudaReductionArrayUSize*sizeof(Real));
    cudaMalloc(&cudaReductionArrayV1, cudaReductionArrayVSize*sizeof(Real));
    cudaMalloc(&cudaReductionArrayV2, cudaReductionArrayVSize*sizeof(Real));
    cudaMalloc(&cudaReductionArrayW1, cudaReductionArrayWSize*sizeof(Real));
    cudaMalloc(&cudaReductionArrayW2, cudaReductionArrayWSize*sizeof(Real));
    cudaMalloc(&cudaReductionArrayResidual1, cudaReductionArrayResidualSize1*sizeof(Real));
    cudaMalloc(&cudaReductionArrayResidual2, cudaReductionArrayResidualSize2*sizeof(Real));
    cudaMalloc(&cudaReductionArrayNumCells1, cudaReductionArrayResidualSize1*sizeof(unsigned int));
    cudaMalloc(&cudaReductionArrayNumCells2, cudaReductionArrayResidualSize2*sizeof(unsigned int));

    // Copy the content of U, V, W, P, T and Flag to the internal representation.
    cudaMemcpy(this->U, U, sizeof(Real)*(imax+1)*(jmax+2)*(kmax+2), cudaMemcpyHostToDevice);
    cudaMemcpy(this->V, V, sizeof(Real)*(imax+2)*(jmax+1)*(kmax+2), cudaMemcpyHostToDevice);
    cudaMemcpy(this->W, W, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+1), cudaMemcpyHostToDevice);
    cudaMemcpy(this->P, P, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), cudaMemcpyHostToDevice);
    cudaMemcpy(this->T, T, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), cudaMemcpyHostToDevice);
    cudaMemcpy(this->Flag, Flag, sizeof(unsigned int)*(imax+2)*(jmax+2)*(kmax+2), cudaMemcpyHostToDevice);
}

CfdSolverCuda::~CfdSolverCuda() {
    cudaFree(U);
    cudaFree(V);
    cudaFree(W);
    cudaFree(P);
    cudaFree(P_temp);
    cudaFree(T);
    cudaFree(T_temp);
    cudaFree(F);
    cudaFree(G);
    cudaFree(H);
    cudaFree(RS);
    cudaFree(Flag);

    cudaFree(cudaReductionArrayU1);
    cudaFree(cudaReductionArrayU2);
    cudaFree(cudaReductionArrayV1);
    cudaFree(cudaReductionArrayV2);
    cudaFree(cudaReductionArrayW1);
    cudaFree(cudaReductionArrayW2);
    cudaFree(cudaReductionArrayResidual1);
    cudaFree(cudaReductionArrayResidual2);
    cudaFree(cudaReductionArrayNumCells1);
    cudaFree(cudaReductionArrayNumCells2);
}

void CfdSolverCuda::setBoundaryValues() {
    setBoundaryValuesCuda(
            T_h, T_c, imax, jmax, kmax, blockSizeX, blockSizeY, blockSizeZ, U, V, W, T, Flag);
}

void CfdSolverCuda::setBoundaryValuesScenarioSpecific() {
    setBoundaryValuesScenarioSpecificCuda(
            scenarioName, imax, jmax, kmax, blockSizeX, blockSizeY, blockSizeZ, U, V, W, Flag);
}

Real CfdSolverCuda::calculateDt() {
    calculateDtCuda(
            Re, Pr, tau, dt, dx, dy, dz, imax, jmax, kmax, blockSize1D, U, V, W,
            cudaReductionArrayU1, cudaReductionArrayU2,
            cudaReductionArrayV1, cudaReductionArrayV2,
            cudaReductionArrayW1, cudaReductionArrayW2,
            useTemperature);
    return dt;
}


void CfdSolverCuda::calculateTemperature() {
    Real *temp = T;
    T = T_temp;
    T_temp = temp;
    dim3 dimBlock(blockSizeX, blockSizeY, blockSizeZ);
    dim3 dimGrid(iceil(kmax,dimBlock.x),iceil(jmax,dimBlock.y),iceil(imax,dimBlock.z));
    calculateTemperatureCuda<<<dimGrid,dimBlock>>>(
            Re, Pr, alpha, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, T, T_temp, Flag);
}

void CfdSolverCuda::calculateFgh() {
    calculateFghCuda(
            Re, GX, GY, GZ, alpha, beta, dt, dx, dy, dz, imax, jmax, kmax,
            blockSizeX, blockSizeY, blockSizeZ, U, V, W, T, F, G, H, Flag);
}

void CfdSolverCuda::calculateRs() {
    dim3 dimBlock(blockSizeX, blockSizeY, blockSizeZ);
    dim3 dimGrid(iceil(kmax,dimBlock.x),iceil(jmax,dimBlock.y),iceil(imax,dimBlock.z));
    calculateRsCuda<<<dimGrid,dimBlock>>>(dt, dx, dy, dz, imax, jmax, kmax, F, G, H, RS);
}


void CfdSolverCuda::executeSorSolver() {
    sorSolverCuda(
            omg, eps, itermax, linearSystemSolverType, shallWriteOutput,
            dx, dy, dz, imax, jmax, kmax,
            blockSizeX, blockSizeY, blockSizeZ, blockSize1D, P, P_temp, RS, Flag,
            cudaReductionArrayResidual1, cudaReductionArrayResidual2,
            cudaReductionArrayNumCells1, cudaReductionArrayNumCells2);
}

void CfdSolverCuda::calculateUvw() {
    dim3 dimBlock(blockSizeX, blockSizeY, blockSizeZ);
    dim3 dimGrid(iceil(kmax,dimBlock.x),iceil(jmax,dimBlock.y),iceil(imax,dimBlock.z));
    calculateUvwCuda<<<dimGrid,dimBlock>>>(dt, dx, dy, dz, imax, jmax, kmax, U, V, W, F, G, H, P, Flag);
}

void CfdSolverCuda::getDataForOutput(Real *U, Real *V, Real *W, Real *P, Real *T) {
    // Copy the content of U, V, W, P, T in the internal representation to the specified output arrays.
    cudaMemcpy(U, this->U, sizeof(Real)*(imax+1)*(jmax+2)*(kmax+2), cudaMemcpyDeviceToHost);
    cudaMemcpy(V, this->V, sizeof(Real)*(imax+2)*(jmax+1)*(kmax+2), cudaMemcpyDeviceToHost);
    cudaMemcpy(W, this->W, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+1), cudaMemcpyDeviceToHost);
    cudaMemcpy(P, this->P, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), cudaMemcpyDeviceToHost);
    cudaMemcpy(T, this->T, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), cudaMemcpyDeviceToHost);
}
