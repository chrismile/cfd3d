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

__global__ void set_x_y_planes_pressure_boundaries(
    int imax, int jmax, int kmax, Real *P){

        int i = blockIdx.x + threadIdx.y + 1;
        int j = blockIdx.y + threadIdx.x + 1;
        
        // Set the boundary values for the pressure on the x-y-planes.
        if (i <= imax && j <= jmax){
                P[IDXP(i,j,0)] = P[IDXP(i,j,1)];
                P[IDXP(i,j,kmax+1)] = P[IDXP(i,j,kmax)];
            }
    }

__global__ void set_x_z_planes_pressure_boundaries(
    int imax, int jmax, int kmax, Real *P){

        int i = blockIdx.x + threadIdx.y + 1;
        int k = blockIdx.y + threadIdx.x + 1;
        // Set the boundary values for the pressure on the x-z-planes.
        if (i <= imax && k<= jmax){
            P[IDXP(i,0,k)] = P[IDXP(i,1,k)];
            P[IDXP(i,jmax+1,k)] = P[IDXP(i,jmax,k)];
        }
    }

__global__ void set_y_z_planes_pressure_boundaries(
    int imax, int jmax, int kmax, Real *P){

        int j = blockIdx.x + threadIdx.y + 1;
        int k = blockIdx.y + threadIdx.x + 1;
        // Set the boundary values for the pressure on the y-z-planes.
        if (j <= imax && k<= jmax){
            P[IDXP(0,j,k)] = P[IDXP(1,j,k)];
            P[IDXP(imax+1,j,k)] = P[IDXP(imax,j,k)];
        }
    }

__global__ void sorSolverIterationCuda(
        Real omg, Real dx, Real dy, Real dz, Real coeff, int imax, int jmax, int kmax,
        Real *P, Real *P_temp, Real *RS, FlagType *Flag, Real &residual) {

        int i = blockIdx.x + 1;
        int j = blockIdx.y + threadIdx.y + 1;
        int k = blockIdx.z + threadIdx.x + 1;   
    
        // Now start with the actual SOR iteration.
    #ifdef SOR_GAUSS_SEIDL
        if (i <= imax && j <= jmax && k <= kmax){
            if (isFluid(Flag[IDXFLAG(i,j,k)])){
                P[IDXP(i,j,k)] = (Real(1.0) - omg)*P[IDXP(i,j,k)] + coeff *
                        ((P[IDXP(i+1,j,k)]+P[IDXP(i-1,j,k)])/(dx*dx)
                        + (P[IDXP(i,j+1,k)]+P[IDXP(i,j-1,k)])/(dy*dy)
                        + (P[IDXP(i,j,k+1)]+P[IDXP(i,j,k-1)])/(dz*dz)
                        - RS[IDXRS(i,j,k)]);
            }
        }
    #endif
    #ifdef SOR_JACOBI
        //#pragma omp parallel for
        if (i <= imax && j <= jmax && k <= kmax){
            P[IDXP(i,j,k)] = (Real(1.0) - omg)*P_temp[IDXP(i,j,k)] + coeff *
                    ((P_temp[IDXP(i+1,j,k)]+P_temp[IDXP(i-1,j,k)])/(dx*dx)
                        + (P_temp[IDXP(i,j+1,k)]+P_temp[IDXP(i,j-1,k)])/(dy*dy)
                        + (P_temp[IDXP(i,j,k+1)]+P_temp[IDXP(i,j,k-1)])/(dz*dz)
                        - RS[IDXRS(i,j,k)]);
        }
    #endif
    #ifdef SOR_HYBRID
        #pragma omp parallel for
        if (i <= imax && j <= jmax && k <= kmax){
            // Just use Jacobi scheme in i direction, as we have only parallelized the outer loop.
            P[IDXP(i,j,k)] = (Real(1.0) - omg)*P[IDXP(i,j,k)] + coeff *
                    ((P_temp[IDXP(i+1,j,k)]+P_temp[IDXP(i-1,j,k)])/(dx*dx)
                        + (P[IDXP(i,j+1,k)]+P[IDXP(i,j-1,k)])/(dy*dy)
                        + (P[IDXP(i,j,k+1)]+P[IDXP(i,j,k-1)])/(dz*dz)
                        - RS[IDXRS(i,j,k)]);
        }
    #endif
    
        // Compute the residual.
        // residual = 0;
        // #pragma omp parallel for reduction(+: residual)
        // for (int i = 1; i <= imax; i++) {
        //     for (int j = 1; j <= jmax; j++) {
        //         for (int k = 1; k <= kmax; k++) {
        //             if (isFluid(Flag[IDXFLAG(i,j,k)])){
        //                 residual += SQR(
        //                            (P[IDXP(i+1,j,k)] - Real(2.0)*P[IDXP(i,j,k)] + P[IDXP(i-1,j,k)])/(dx*dx)
        //                          + (P[IDXP(i,j+1,k)] - Real(2.0)*P[IDXP(i,j,k)] + P[IDXP(i,j-1,k)])/(dy*dy)
        //                          + (P[IDXP(i,j,k+1)] - Real(2.0)*P[IDXP(i,j,k)] + P[IDXP(i,j,k-1)])/(dz*dz)
        //                          - RS[IDXRS(i,j,k)]
        //                 );
        //             }
        //         }
        //     }
        // }
        // TODO
    
        // The residual is normalized by dividing by the total number of fluid cells.
        //residual = std::sqrt(residual/(imax*jmax*kmax));
}

void sorSolverCuda(Real omg, Real eps, int itermax,
        Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *P, Real *P_temp, Real *RS, FlagType *Flag) {
    const Real coeff = omg / (2.0 * (1.0 / (dx*dx) + 1.0 / (dy*dy) + 1.0 / (dz*dz)));
    Real residual = 1e9;
    int it = 0;
    while (it < itermax && residual > eps) {
        dim3 dimBlock(32,32);
        dim3 dimGrid_x_y(iceil(imax,dimBlock.y),iceil(jmax,dimBlock.x));
        set_x_y_planes_pressure_boundaries<<<dimGrid_x_y,dimBlock>>>(imax, jmax, kmax, P);

        dim3 dimGrid_x_z(iceil(imax,dimBlock.y),iceil(kmax,dimBlock.x));
        set_x_z_planes_pressure_boundaries<<<dimGrid_x_z,dimBlock>>>(imax, jmax, kmax, P);

        dim3 dimGrid_y_z(iceil(jmax,dimBlock.y),iceil(kmax,dimBlock.x));
        set_y_z_planes_pressure_boundaries<<<dimGrid_y_z,dimBlock>>>(imax, jmax, kmax, P);

#if defined(SOR_JACOBI) || defined(SOR_HYBRID)
        cudaMemcpy(P_temp, P, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), cudaMemcpyDeviceToDevice);
#endif

        dim3 dimGrid(iceil(imax,dimBlock.z),iceil(jmax,dimBlock.y),iceil(kmax,dimBlock.x));
        sorSolverIterationCuda<<<dimGrid,dimBlock>>>(omg, dx, dy, dz, coeff, imax, jmax, kmax, P, P_temp, RS, Flag, residual);

        it++;
    }

    if (residual > eps && it == itermax || std::isnan(residual)) {
        std::cerr << "\nSOR solver reached maximum number of iterations without converging." << std::endl;
    }
}
