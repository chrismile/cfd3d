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
#include <cmath>
#include <cstring>
#include "../Flag.hpp"
#include "SorSolverCpp.hpp"

void sorSolverIterationCpp(
        Real omg, Real dx, Real dy, Real dz, Real coeff, int imax, int jmax, int kmax,
        LinearSystemSolverType linearSystemSolverType,
        Real *P, Real *P_temp, Real *RS, FlagType *Flag, Real &residual) {
    // Set the boundary values for the pressure on the x-y-planes.
    #pragma omp parallel for
    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            P[IDXP(i,j,0)] = P[IDXP(i,j,1)];
            P[IDXP(i,j,kmax+1)] = P[IDXP(i,j,kmax)];
        }
    }

    // Set the boundary values for the pressure on the x-z-planes.
    #pragma omp parallel for
    for (int i = 1; i <= imax; i++) {
        for (int k = 1; k <= kmax; k++) {
            P[IDXP(i,0,k)] = P[IDXP(i,1,k)];
            P[IDXP(i,jmax+1,k)] = P[IDXP(i,jmax,k)];
        }
    }

    // Set the boundary values for the pressure on the y-z-planes.
    #pragma omp parallel for
    for (int j = 1; j <= jmax; j++) {
        for (int k = 1; k <= kmax; k++) {
            P[IDXP(0,j,k)] = P[IDXP(1,j,k)];
            P[IDXP(imax+1,j,k)] = P[IDXP(imax,j,k)];
        }
    }

    // Boundary values for arbitrary geometries.
    #pragma omp parallel for
    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax; k++) {
                int numDirectFlag = 0;
                Real P_temp = Real(0);

                if (!isFluid(Flag[IDXFLAG(i, j, k)])) {
                    if (B_R(Flag[IDXFLAG(i, j, k)])) {
                        P_temp += P[IDXP(i + 1, j, k)];
                        numDirectFlag++;
                    }

                    if (B_L(Flag[IDXFLAG(i, j, k)])) {
                        P_temp += P[IDXP(i - 1, j, k)];
                        numDirectFlag++;
                    }

                    if (B_U(Flag[IDXFLAG(i, j, k)])) {
                        P_temp += P[IDXP(i, j + 1, k)];
                        numDirectFlag++;
                    }

                    if (B_D(Flag[IDXFLAG(i, j, k)])) {
                        P_temp += P[IDXP(i, j - 1, k)];
                        numDirectFlag++;
                    }

                    if (B_B(Flag[IDXFLAG(i, j, k)])) {
                        P_temp += P[IDXP(i, j, k - 1)];
                        numDirectFlag++;
                    }

                    if (B_F(Flag[IDXFLAG(i, j, k)])) {
                        P_temp += P[IDXP(i, j, k + 1)];
                        numDirectFlag++;
                    }

                    P[IDXP(i, j, k)] = P_temp / Real(numDirectFlag);
                }
            }     
        }
    }



    if (linearSystemSolverType == LINEAR_SOLVER_JACOBI) {
        // Create a copy of the current state of the pressure array.
        // A parallel loop is potentially faster than 'memcpy' for large domains.
        #pragma omp parallel for
        for (int i = 0; i <= imax+1; i++) {
            for (int j = 0; j <= jmax+1; j++) {
                for (int k = 0; k <= kmax+1; k++) {
                    P_temp[IDXP(i, j, k)] = P[IDXP(i, j, k)];
                }
            }
        }
    }


    if (linearSystemSolverType == LINEAR_SOLVER_SOR || linearSystemSolverType == LINEAR_SOLVER_GAUSS_SEIDEL) {
        for (int i = 1; i <= imax; i++) {
            for (int j = 1; j <= jmax; j++) {
                for (int k = 1; k <= kmax; k++) {
                    if (isFluid(Flag[IDXFLAG(i,j,k)])){
                        P[IDXP(i,j,k)] = (Real(1.0) - omg)*P[IDXP(i,j,k)] + coeff *
                                ((P[IDXP(i+1,j,k)]+P[IDXP(i-1,j,k)])/(dx*dx)
                                + (P[IDXP(i,j+1,k)]+P[IDXP(i,j-1,k)])/(dy*dy)
                                + (P[IDXP(i,j,k+1)]+P[IDXP(i,j,k-1)])/(dz*dz)
                                - RS[IDXRS(i,j,k)]);
                    }
                }
            }
        }
    } else if (linearSystemSolverType == LINEAR_SOLVER_SOR_PARALLEL
            || linearSystemSolverType == LINEAR_SOLVER_GAUSS_SEIDEL_PARALLEL) {
        #pragma omp parallel for
        for (int i = 1; i <= imax; i++) {
            for (int j = 1; j <= jmax; j++) {
                for (int k = 1; k <= kmax; k++) {
                    if (isFluid(Flag[IDXFLAG(i,j,k)])){
                        P_temp[IDXP(i,j,k)] = (Real(1.0) - omg)*P[IDXP(i,j,k)] + coeff *
                                ((P[IDXP(i+1,j,k)])/(dx*dx)
                                + (P[IDXP(i,j+1,k)])/(dy*dy)
                                + (P[IDXP(i,j,k+1)])/(dz*dz)
                                - RS[IDXRS(i,j,k)]);
                    }
                }
            }
        }

        for (int i = 1; i <= imax; i++) {
            for (int j = 1; j <= jmax; j++) {
                for (int k = 1; k <= kmax; k++) {
                    if (isFluid(Flag[IDXFLAG(i,j,k)])){
                        P[IDXP(i,j,k)] = P_temp[IDXP(i,j,k)] + coeff *
                                ((P[IDXP(i-1,j,k)])/(dx*dx)
                                 + (P[IDXP(i,j-1,k)])/(dy*dy)
                                 + (P[IDXP(i,j,k-1)])/(dz*dz));
                    }
                }
            }
        }
    } else if (linearSystemSolverType == LINEAR_SOLVER_JACOBI) {
        #pragma omp parallel for
        for (int i = 1; i <= imax; i++) {
            for (int j = 1; j <= jmax; j++) {
                for (int k = 1; k <= kmax; k++) {
                    if (isFluid(Flag[IDXFLAG(i,j,k)])) {
                        P[IDXP(i,j,k)] = (Real(1.0) - omg)*P_temp[IDXP(i,j,k)] + coeff *
                                ((P_temp[IDXP(i+1,j,k)]+P_temp[IDXP(i-1,j,k)])/(dx*dx)
                                 + (P_temp[IDXP(i,j+1,k)]+P_temp[IDXP(i,j-1,k)])/(dy*dy)
                                 + (P_temp[IDXP(i,j,k+1)]+P_temp[IDXP(i,j,k-1)])/(dz*dz)
                                 - RS[IDXRS(i,j,k)]);
                    }
                }
            }
        }
    }


    // Compute the residual.
    residual = Real(0.0);
    int numFluidCells = 0;
    #pragma omp parallel for reduction(+: residual) reduction(+: numFluidCells)
    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax; k++) {
                if (isFluid(Flag[IDXFLAG(i,j,k)])){
                    residual += SQR(
                               (P[IDXP(i+1,j,k)] - Real(2.0)*P[IDXP(i,j,k)] + P[IDXP(i-1,j,k)])/(dx*dx)
                             + (P[IDXP(i,j+1,k)] - Real(2.0)*P[IDXP(i,j,k)] + P[IDXP(i,j-1,k)])/(dy*dy)
                             + (P[IDXP(i,j,k+1)] - Real(2.0)*P[IDXP(i,j,k)] + P[IDXP(i,j,k-1)])/(dz*dz)
                             - RS[IDXRS(i,j,k)]
                    );
                    numFluidCells++;
                }
            }
        }
    }

    // The residual is normalized by dividing by the total number of fluid cells.
    residual = std::sqrt(residual/numFluidCells);
}

void sorSolverCpp(
        Real omg, Real eps, int itermax, LinearSystemSolverType linearSystemSolverType,
        Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *P, Real *P_temp, Real *RS, FlagType *Flag) {
    if (linearSystemSolverType == LINEAR_SOLVER_SOR || linearSystemSolverType == LINEAR_SOLVER_SOR_PARALLEL) {
        // Successive over-relaxation based on Gauss-Seidl. A factor of 1.5 proved to give the best results here.
        omg = 1.5;
    } else {
        // A method named JOR (Jacobi over-relaxation) with omega != 1 exists, but doesn't converge for this problem.
        omg = 1.0;
    }

    const Real coeff = omg / (Real(2.0) * (Real(1.0) / (dx*dx) + Real(1.0) / (dy*dy) + Real(1.0) / (dz*dz)));
    Real residual = Real(1e9);
    int it = 0;

    while (it < itermax && residual > eps) {
        sorSolverIterationCpp(
                omg, dx, dy, dz, coeff, imax, jmax, kmax, linearSystemSolverType,
                P, P_temp, RS, Flag, residual);
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
