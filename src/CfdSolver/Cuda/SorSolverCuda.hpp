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

#ifndef CFD3D_SORSOLVERCUDA_HPP
#define CFD3D_SORSOLVERCUDA_HPP

#include "Defines.hpp"

/*
* Sets the x_y plane boundary values for the pressure
*/
__global__ void set_x_y_planes_pressure_boundaries(
    int imax, int jmax, int kmax, Real *P);

__global__ void set_x_z_planes_pressure_boundaries(
    int imax, int jmax, int kmax, Real *P);

__global__ void set_y_z_planes_pressure_boundaries(
    int imax, int jmax, int kmax, Real *P);

__global__ void copy_pressure(
    int imax, int jmax, int kmax,Real *P, Real *P_temp);

/**
 * Uses an SOR solver to compute the updated pressure values using the pressure poisson equation (PPE).
 */
void sorSolverCuda(
        Real omg, Real eps, int itermax, LinearSystemSolverType linearSystemSolverType,
        Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        int blockSizeX, int blockSizeY, int blockSizeZ, int blockSize1D,
        Real *P, Real *P_temp, Real *RS, FlagType *Flag,
        Real *cudaReductionArrayResidual1, Real *cudaReductionArrayResidual2,
        unsigned int *cudaReductionArrayNumCells1, unsigned int *cudaReductionArrayNumCells2);

#endif //CFD3D_SORSOLVERCUDA_HPP
