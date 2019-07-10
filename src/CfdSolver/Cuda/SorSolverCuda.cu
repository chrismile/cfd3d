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

void sorSolverIterationCuda(
        Real omg, Real dx, Real dy, Real dz, Real coeff, int imax, int jmax, int kmax,
        Real *P, Real *P_temp, Real *RS, FlagType *Flag, Real &residual) {
    // TODO
}

void sorSolverCuda(
        Real omg, Real eps, int itermax, LinearSystemSolverType linearSystemSolverType,
        Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        int blockSizeX, int blockSizeY, int blockSizeZ, int blockSize1D,
        Real *P, Real *P_temp, Real *RS, FlagType *Flag) {
    if (linearSystemSolverType == LINEAR_SOLVER_SOR || linearSystemSolverType == LINEAR_SOLVER_SOR_PARALLEL) {
        // Successive over-relaxation based on Gauss-Seidl. A factor of 1.5 proved to give the best results here.
        omg = 1.5;
    } else {
        // A method named JOR (Jacobi over-relaxation) with omega != 1 exists, but doesn't converge for this problem.
        omg = 1.0;
    }

    const Real coeff = omg / (2.0 * (1.0 / (dx*dx) + 1.0 / (dy*dy) + 1.0 / (dz*dz)));
    Real residual = 1e9;
    int it = 0;
    while (it < itermax && residual > eps) {
        cudaMemcpy(P_temp, P, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), cudaMemcpyDeviceToDevice);
        sorSolverIterationCuda(omg, dx, dy, dz, coeff, imax, jmax, kmax, P, P_temp, RS, Flag, residual);
        it++;
    }

    if (residual > eps && it == itermax) {
        std::cerr << "\nSOR solver reached maximum number of iterations without converging." << std::endl;
    }
}
