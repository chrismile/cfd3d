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
#include "BoundaryValuesCuda.hpp"
#include "UvwCuda.hpp"
#include "SorSolverCuda.hpp"
#include "CfdSolverCuda.hpp"
#include "../Cpp/SorSolverCpp.hpp"
#include "../Cpp/BoundaryValuesCpp.hpp"
#include "../Cpp/UvwCpp.hpp"


void CfdSolverCuda::initialize(const std::string &scenarioName,
        Real Re, Real Pr, Real omg, Real eps, int itermax, Real alpha, Real beta, Real dt, Real tau,
        Real GX, Real GY, Real GZ, bool useTemperature, Real T_h, Real T_c,
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz,
        Real *U, Real *V, Real *W, Real *P, Real *T, uint32_t *Flag) {
    this->scenarioName = scenarioName;
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
}

void CfdSolverCuda::setBoundaryValues() {
    setBoundaryValuesCpp(T_h, T_c, imax, jmax, kmax, U, V, W, T, Flag);
}

void CfdSolverCuda::setBoundaryValuesScenarioSpecific() {
    setBoundaryValuesScenarioSpecificCpp(scenarioName, imax, jmax, kmax, U, V, W, Flag);
}

Real CfdSolverCuda::calculateDt() {
    calculateDtCpp(Re, Pr, tau, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, useTemperature);
    return dt;
}


void CfdSolverCuda::calculateTemperature() {
    Real *temp = T;
    T = T_temp;
    T_temp = temp;
    calculateTemperatureCpp(Re, Pr, alpha, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, T, T_temp, Flag);
}

void CfdSolverCuda::calculateFgh() {
    calculateFghCuda<<<dimGrid,dimBlock>>>(Re, GX, GY, GZ, alpha, beta, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, T, F, G, H, Flag)
}

void CfdSolverCuda::calculateRs() {
    calculateRsCpp<<<dimGrid,dimBlock>>>(dt, dx, dy, dz, imax, jmax, kmax, F, G, H, RS);
}


void CfdSolverCuda::executeSorSolver() {
    sorSolverCpp(omg, eps, itermax, dx, dy, dz, imax, jmax, kmax, P, P_temp, RS, Flag);
}

void CfdSolverCuda::calculateUvw() {
    calculateUvwCpp(dt, dx, dy, dz, imax, jmax, kmax, U, V, W, F, G, H, P, Flag);
}

void CfdSolverCuda::getDataForOutput(Real *U, Real *V, Real *W, Real *P, Real *T) {
    // Copy the content of U, V, W, P, T in the internal representation to the specified output arrays.
    cudaMemcpy(U, this->U, sizeof(Real)*(imax+1)*(jmax+2)*(kmax+2), cudaMemcpyDeviceToHost);
    cudaMemcpy(V, this->V, sizeof(Real)*(imax+2)*(jmax+1)*(kmax+2), cudaMemcpyDeviceToHost);
    cudaMemcpy(W, this->W, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+1), cudaMemcpyDeviceToHost);
    cudaMemcpy(P, this->P, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), cudaMemcpyDeviceToHost);
    cudaMemcpy(T, this->T, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), cudaMemcpyDeviceToHost);
    cudaMemcpy(Flag, this->Flag, sizeof(unsigned int)*(imax+2)*(jmax+2)*(kmax+2), cudaMemcpyDeviceToHost);
}
