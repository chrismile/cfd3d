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
#include "MpiHelpers.hpp"
#include "DefinesMpi.hpp"
#include "SorSolverMpi.hpp"

// Three possible modes: Gauss-Seidl, Jacobi, Gauss-Seidl/Jacobi hybrid
#define SOR_GAUSS_SEIDL
//#define SOR_GAUSS_SEIDL_PARALLEL
//#define SOR_JACOBI
//#define SOR_HYBRID

void sorSolverIterationMpi(
        Real omg, Real dx, Real dy, Real dz, Real coeff, int imax, int jmax, int kmax,
        int il, int iu, int jl, int ju, int kl, int ku,
        int rankL, int rankR, int rankD, int rankU, int rankB, int rankF, Real *bufSend, Real *bufRecv,
        Real *P, Real *P_temp, Real *RS, FlagType *Flag, Real &residual) {
    // Set the boundary values for the pressure on the x-y-planes.
    if (kl == 1) {
        for (int i = 1; i <= imax; i++) {
            for (int j = 1; j <= jmax; j++) {
                P[IDXP(i, j, 0)] = P[IDXP(i, j, 1)];
            }
        }
    }
    if (ku == kmax) {
        for (int i = 1; i <= imax; i++) {
            for (int j = 1; j <= jmax; j++) {
                P[IDXP(i, j, kmax + 1)] = P[IDXP(i, j, kmax)];
            }
        }
    }

    // Set the boundary values for the pressure on the x-z-planes.
    if (jl == 1) {
        for (int i = 1; i <= imax; i++) {
            for (int k = 1; k <= kmax; k++) {
                P[IDXP(i,0,k)] = P[IDXP(i,1,k)];
            }
        }
    }
    if (ju == jmax) {
        for (int i = 1; i <= imax; i++) {
            for (int k = 1; k <= kmax; k++) {
                P[IDXP(i,jmax+1,k)] = P[IDXP(i,jmax,k)];
            }
        }
    }

    // Set the boundary values for the pressure on the y-z-planes.
    if (il == 1) {
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax; k++) {
                P[IDXP(0,j,k)] = P[IDXP(1,j,k)];
            }
        }
    }
    if (iu == imax) {
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax; k++) {
                P[IDXP(imax+1,j,k)] = P[IDXP(imax,j,k)];
            }
        }
    }

    for (int i = il; i <= iu; i++) {
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
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



#if defined(SOR_JACOBI) || defined(SOR_HYBRID)
    //memcpy(P_temp, P, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2));
    // Using multiple threads for the copy is faster for large amounts of data.
    for (int i = il-1; i <= iu+1; i++) {
        for (int j = jl-1; j <= ju+1; j++) {
            for (int k = kl-1; k <= ku+1; k++) {
                P_temp[IDXP(i, j, k)] = P[IDXP(i, j, k)];
            }
        }
    }
#endif


    // Now start with the actual SOR iteration.
#ifdef SOR_GAUSS_SEIDL
    for (int i = il; i <= iu; i++) {
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
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
#endif
#ifdef SOR_GAUSS_SEIDL_PARALLEL
    for (int i = il; i <= iu; i++) {
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
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

    for (int i = il; i <= iu; i++) {
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                if (isFluid(Flag[IDXFLAG(i,j,k)])){
                    P[IDXP(i,j,k)] = P_temp[IDXP(i,j,k)] + coeff *
                            ((P[IDXP(i-1,j,k)])/(dx*dx)
                            + (P[IDXP(i,j-1,k)])/(dy*dy)
                            + (P[IDXP(i,j,k-1)])/(dz*dz));
                }
            }
        }
    }
#endif
#ifdef SOR_JACOBI
    for (int i = il; i <= iu; i++) {
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
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
#endif
#ifdef SOR_HYBRID
    #pragma omp parallel for
    for (int i = il; i <= iu; i++) {
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                // Just use Jacobi scheme in i direction, as we have only parallelized the outer loop.
                P[IDXP(i,j,k)] = (Real(1.0) - omg)*P[IDXP(i,j,k)] + coeff *
                        ((P_temp[IDXP(i+1,j,k)]+P_temp[IDXP(i-1,j,k)])/(dx*dx)
                         + (P[IDXP(i,j+1,k)]+P[IDXP(i,j-1,k)])/(dy*dy)
                         + (P[IDXP(i,j,k+1)]+P[IDXP(i,j,k-1)])/(dz*dz)
                         - RS[IDXRS(i,j,k)]);
            }
        }
    }
#endif

    MPI_Status *status = 0;
    mpiExchangeCellData(P, il, iu, jl, ju, kl, ku, rankL, rankR, rankD, rankU, rankB, rankF, bufSend, bufRecv, status);

    // Compute the residual.
    residual = 0;
    int numFluidCells = 0;
    #pragma omp parallel for reduction(+: residual) reduction(+: numFluidCells)
    for (int i = il; i <= iu; i++) {
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
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
    MPI_Allreduce(MPI_IN_PLACE, &numFluidCells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &residual, 1, MPI_REAL_CFD3D, MPI_SUM, MPI_COMM_WORLD);
    residual = std::sqrt(residual/numFluidCells);
}

void sorSolverMpi(
        Real omg, Real eps, int itermax,
        Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        int il, int iu, int jl, int ju, int kl, int ku,
        int rankL, int rankR, int rankD, int rankU, int rankB, int rankF, Real *bufSend, Real *bufRecv,
        Real *P, Real *P_temp, Real *RS, FlagType *Flag) {
#if defined(SOR_JACOBI) || defined(SOR_HYBRID)
    omg = 1.0;
#else
    omg = 1.5;
#endif

    const Real coeff = omg / (Real(2.0) * (Real(1.0) / (dx*dx) + Real(1.0) / (dy*dy) + Real(1.0) / (dz*dz)));
    Real residual = Real(1e9);
    int it = 0;

    while (it < itermax && residual > eps) {
        sorSolverIterationMpi(
                omg, dx, dy, dz, coeff, imax, jmax, kmax, il, iu, jl, ju, kl, ku,
                rankL, rankR, rankD, rankU, rankB, rankF, bufSend, bufRecv,
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
