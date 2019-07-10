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

#include "BoundaryValuesCuda.hpp"
#include "../../Defines.hpp"
#include "CudaDefines.hpp"

__global__ void setLeftRightBoundariesCuda(
    Real T_h, Real T_c,
    int imax, int jmax, int kmax,
    Real *U, Real *V, Real *W, Real *T,
    FlagType *Flag) {

    int j = blockIdx.y + threadIdx.y + 1;
    int k = blockIdx.x + threadIdx.x + 1;
    // Set the boundary values for the pressure on the y-z-planes.
    if (j <= imax && k<= jmax){
        // Left wall
        if (isNoSlip(Flag[IDXFLAG(0,j,k)])) {
            U[IDXU(0,j,k)] = 0.0;
            V[IDXV(0,j,k)] = -V[IDXV(1,j,k)];
            W[IDXW(0,j,k)] = -W[IDXW(1,j,k)];
        }
        else if (isFreeSlip(Flag[IDXFLAG(0,j,k)])) {
            U[IDXU(0,j,k)] = 0.0;
            V[IDXV(0,j,k)] = V[IDXV(1,j,k)];
            W[IDXW(0,j,k)] = W[IDXW(1,j,k)];
        }
        else if (isOutflow(Flag[IDXFLAG(0,j,k)])) {
            U[IDXU(0,j,k)] = U[IDXU(1,j,k)];
            V[IDXV(0,j,k)] = V[IDXV(1,j,k)];
            W[IDXW(0,j,k)] = W[IDXW(1,j,k)];
        }
        // Right wall
        if (isNoSlip(Flag[IDXFLAG(imax+1,j,k)])) {
            U[IDXU(imax,j,k)] = 0.0;
            V[IDXV(imax+1,j,k)] = -V[IDXV(imax,j,k)];
            W[IDXW(imax+1,j,k)] = -W[IDXW(imax,j,k)];
        }
        else if (isFreeSlip(Flag[IDXFLAG(imax+1,j,k)])) {
            U[IDXU(imax,j,k)] = 0.0;
            V[IDXV(imax+1,j,k)] = V[IDXV(imax,j,k)];
            W[IDXW(imax+1,j,k)] = W[IDXW(imax,j,k)];
        }
        else if (isOutflow(Flag[IDXFLAG(imax+1,j,k)])) {
            U[IDXU(imax,j,k)] = U[IDXU(imax-1,j,k)];
            V[IDXV(imax+1,j,k)] = V[IDXV(imax,j,k)];
            W[IDXW(imax+1,j,k)] = W[IDXW(imax,j,k)];
        }

        // Left boundary T
        if (isHot(Flag[IDXFLAG(0,j,k)])) {
            T[IDXT(0,j,k)] = 2 * T_h - T[IDXT(1,j,k)];
        } else if (isCold(Flag[IDXFLAG(0,j,k)])) {
            T[IDXT(0,j,k)] = 2 * T_c - T[IDXT(1,j,k)];
        } else {
            T[IDXT(0,j,k)] = T[IDXT(1,j,k)];
        }
        
        // Right boundary T
        if (isHot(Flag[IDXFLAG(imax+1,j,k)])) {
            T[IDXT(imax+1,j,k)] = 2 * T_h - T[IDXT(imax,j,k)];
        }  else if (isCold(Flag[IDXFLAG(imax+1,j,k)])) {
            T[IDXT(imax+1,j,k)] = 2 * T_c - T[IDXT(imax,j,k)];
        } else {
            T[IDXT(imax+1,j,k)] = T[IDXT(imax,j,k)];
        }            
    }
}

__global__ void setDownUpBoundariesCuda(
    Real T_h, Real T_c,
    int imax, int jmax, int kmax,
    Real *U, Real *V, Real *W, Real *T,
    FlagType *Flag) {
    int i = blockIdx.y + threadIdx.y + 1;
    int k = blockIdx.x + threadIdx.x + 1;
    // Set the boundary values for the pressure on the x-z-planes.
    if (i <= imax && k<= jmax){
        // Down wall
        if (isNoSlip(Flag[IDXFLAG(i,0,k)])) {
            U[IDXU(i,0,k)] = -U[IDXU(i,1,k)];
            V[IDXV(i,0,k)] = 0.0;
            W[IDXW(i,0,k)] = -W[IDXW(i,1,k)];
        } else if (isFreeSlip(Flag[IDXFLAG(i,0,k)])) {
            U[IDXU(i,0,k)] = U[IDXU(i,1,k)];
            V[IDXV(i,0,k)] = 0.0;
            W[IDXW(i,0,k)] = W[IDXW(i,1,k)];
        } else if (isOutflow(Flag[IDXFLAG(i,0,k)])) {
            U[IDXU(i,0,k)] = U[IDXU(i,1,k)];
            V[IDXV(i,0,k)] = V[IDXV(i,1,k)];
            W[IDXW(i,0,k)] = W[IDXW(i,1,k)];
        }
        // Up wall
        if (isNoSlip(Flag[IDXFLAG(i,jmax+1,k)])) {
            U[IDXU(i,jmax+1,k)] = -U[IDXU(i,jmax,k)];
            V[IDXV(i,jmax,k)] = 0.0;
            W[IDXW(i,jmax+1,k)] = -W[IDXW(i,jmax,k)];
        } else if (isFreeSlip(Flag[IDXFLAG(i,jmax+1,k)])) {
            U[IDXU(i,jmax+1,k)] = U[IDXU(i,jmax,k)];
            V[IDXV(i,jmax,k)] = 0.0;
            W[IDXW(i,jmax+1,k)] = W[IDXW(i,jmax,k)];
        } else if (isOutflow(Flag[IDXFLAG(i,jmax+1,k)])) {
            U[IDXU(i,jmax+1,k)] = U[IDXU(i,jmax,k)];
            V[IDXV(i,jmax,k)] = V[IDXV(i,jmax-1,k)];
            W[IDXW(i,jmax+1,k)] = W[IDXW(i,jmax,k)];
        }

        // Down boundary T
        if (isHot(Flag[IDXFLAG(i,0,k)])) {
            T[IDXT(i,0,k)] = 2 * T_h - T[IDXT(i,1,k)];
        } else if (isCold(Flag[IDXFLAG(i,0,k)])) {
            T[IDXT(i,0,k)] = 2 * T_c - T[IDXT(i,1,k)];
        } else {
            T[IDXT(i,0,k)] = T[IDXT(i,1,k)];
        }
        
        // Up boundary T
        if (isHot(Flag[IDXFLAG(i,jmax+1,k)])) {
            T[IDXT(i,jmax+1,k)] = 2 * T_h - T[IDXT(i,jmax,k)];
        } else if (isCold(Flag[IDXFLAG(i,jmax+1,k)])) {
            T[IDXT(i,jmax+1,k)] = 2 * T_c - T[IDXT(i,jmax,k)];
        } else {
            T[IDXT(i,jmax+1,k)] = T[IDXT(i,jmax,k)];
        }            
    }
}

__global__ void setFrontBackBoundariesCuda(
    Real T_h, Real T_c,
    int imax, int jmax, int kmax,
    Real *U, Real *V, Real *W, Real *T,
    FlagType *Flag) {
    int i = blockIdx.y + threadIdx.y + 1;
    int j = blockIdx.x + threadIdx.x + 1;
    
    // Set the boundary values for the pressure on the x-y-planes.
    if (i <= imax && j <= jmax){
        // Front wall
        if (isNoSlip(Flag[IDXFLAG(i,j,0)])) {
            U[IDXU(i,j,0)] = -U[IDXU(i,j,1)];
            V[IDXV(i,j,0)] = -V[IDXV(i,j,1)];
            W[IDXW(i,j,0)] = 0.0;
        }
        else if (isFreeSlip(Flag[IDXFLAG(i,j,0)])) {
            U[IDXU(i,j,0)] = U[IDXU(i,j,1)];
            V[IDXV(i,j,0)] = V[IDXV(i,j,1)];
            W[IDXW(i,j,0)] = 0.0;
        }
        else if (isOutflow(Flag[IDXFLAG(i,j,0)])) {
            U[IDXU(i,j,0)] = U[IDXU(i,j,1)];
            V[IDXV(i,j,0)] = V[IDXV(i,j,1)];
            W[IDXW(i,j,0)] = W[IDXW(i,j,1)];
        }
        // Back wall
        if (isNoSlip(Flag[IDXFLAG(i,j,kmax+1)])) {
            U[IDXU(i,j,kmax+1)] = -U[IDXU(i,j,kmax)];
            V[IDXV(i,j,kmax+1)] = -V[IDXV(i,j,kmax)];
            W[IDXW(i,j,kmax)] = 0.0;
        }
        else if (isFreeSlip(Flag[IDXFLAG(i,j,kmax+1)])) {
            U[IDXU(i,j,kmax+1)] = U[IDXU(i,j,kmax)];
            V[IDXV(i,j,kmax+1)] = V[IDXV(i,j,kmax)];
            W[IDXW(i,j,kmax)] = 0.0;
        }
        else if (isOutflow(Flag[IDXFLAG(i,j,kmax+1)])) {
            U[IDXU(i,j,kmax+1)] = U[IDXU(i,j,kmax)];
            V[IDXV(i,j,kmax+1)] = V[IDXV(i,j,kmax)];
            W[IDXW(i,j,kmax)] = W[IDXW(i,j,kmax-1)];
        }

        // Front boundary T
        if (isHot(Flag[IDXFLAG(i,j,0)])) {
            T[IDXT(i,j,0)] = 2 * T_h - T[IDXT(i,j,1)];
        } 
        else if (isCold(Flag[IDXFLAG(i,j,0)])) {
            T[IDXT(i,j,0)] = 2 * T_c - T[IDXT(i,j,1)];
        }
        else {
            T[IDXT(i,j,0)] = T[IDXT(i,j,1)];
        }
        
        // Back boundary T
        if (isHot(Flag[IDXFLAG(i,j,kmax+1)])) {
            T[IDXT(i,j,kmax+1)] = 2 * T_h - T[IDXT(i,j,kmax)];
        } 
        else if (isCold(Flag[IDXFLAG(i,j,kmax+1)])) {
            T[IDXT(i,j,kmax+1)] = 2 * T_c - T[IDXT(i,j,kmax)];
        }
        else {
            T[IDXT(i,j,kmax+1)] = T[IDXT(i,j,kmax)];
        }            
    }
}

void setBoundaryValuesCuda(
        Real T_h, Real T_c,
        int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *T,
        FlagType *Flag) {
    dim3 dimBlock(blockSize,blockSize);
    dim3 dimGrid_x_y(iceil(jmax,dimBlock.x),iceil(kmax,dimBlock.x));
    setFrontBackBoundariesCuda<<<dimGrid_x_y,dimBlock>>>(T_h, T_c, imax, jmax, kmax, U, V, W, T, Flag);

    dim3 dimGrid_x_z(iceil(kmax,dimBlock.x),iceil(imax,dimBlock.y));
    setDownUpBoundariesCuda<<<dimGrid_x_z,dimBlock>>>(T_h, T_c, imax, jmax, kmax, U, V, W, T, Flag);

    dim3 dimGrid_y_z(iceil(kmax,dimBlock.x),iceil(jmax,dimBlock.y));
    setLeftRightBoundariesCuda<<<dimGrid_y_z,dimBlock>>>(T_h, T_c, imax, jmax, kmax, U, V, W, T, Flag);
}

__global__ void setDrivenCavityBoundariesCuda(int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W,
        FlagType *Flag){
    
    int i = blockIdx.y + threadIdx.y + 1;
    int k = blockIdx.x + threadIdx.x + 1;

    if (i <= imax && k<= jmax){
        for (int i = 1; i <= imax; i++) {
            for (int k = 1; k <= kmax; k++) {
                // Upper wall
                U[IDXU(i,jmax+1,k)] = 2.0-U[IDXU(i,jmax,k)];
            }
        }
    }
}

void setBoundaryValuesScenarioSpecificCuda(
        const std::string &scenarioName,
        int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W,
        FlagType *Flag) {

    dim3 dimBlock(blockSize,blockSize);
    if (scenarioName == "driven_cavity") {
        dim3 dimGrid_x_z(iceil(kmax,dimBlock.x),iceil(imax,dimBlock.y));
        setDrivenCavityBoundariesCuda<<<dimBlock, dimGrid_x_z>>>(imax, jmax, kmax, U, V, W, Flag);
    }   
}
