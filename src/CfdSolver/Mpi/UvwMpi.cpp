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

#include <cmath>
#include <algorithm>
#include "Defines.hpp"
#include "UvwMpi.hpp"
#include "MpiHelpers.hpp"
#include "../Flag.hpp"
#include "DefinesMpi.hpp"

void calculateFghMpi(
        Real Re, Real GX, Real GY, Real GZ, Real alpha, Real beta,
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        int il, int iu, int jl, int ju, int kl, int ku,
        Real *U, Real *V, Real *W, Real *T, Real *F, Real *G, Real *H, FlagType *Flag) {
    Real d2u_dx2,d2u_dy2,d2u_dz2,
         d2v_dx2,d2v_dy2,d2v_dz2,
         d2w_dx2,d2w_dy2,d2w_dz2;
    
    Real du2_dx,duv_dy,duw_dz,
         duv_dx,dv2_dy,dvw_dz,
         duw_dx,dvw_dy,dw2_dz;
    
    Real Dx = 1/dx, Dy = 1/dy, Dz = 1/dz;

    #pragma omp parallel for private(d2u_dx2, d2u_dy2, d2u_dz2, du2_dx, duv_dy, duw_dz)
    for (int i = il-1; i <= iu; i++) {
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
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
        }
    }

    #pragma omp parallel for private(d2v_dx2, d2v_dy2, d2v_dz2, duv_dx, dv2_dy, dvw_dz)
    for (int i = il; i <= iu; i++) {
        for (int j = jl-1; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
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
        }
    }

    #pragma omp parallel for private(d2w_dx2, d2w_dy2, d2w_dz2, duw_dx, dvw_dy, dw2_dz)
    for (int i = il; i <= iu; i++) {
        for (int j = jl; j <= ju; j++) {
            for (int k = kl-1; k <= ku; k++) {
                if(isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i,j,k+1)])){
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
    }

    if (il == 1) {
        #pragma omp parallel for
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                F[IDXF(0,j,k)] = U[IDXU(0,j,k)];
            }
        }
    }
    if (iu == imax) {
        #pragma omp parallel for
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                F[IDXF(imax,j,k)] = U[IDXU(imax,j,k)];
            }
        }
    }

    if (jl == 1) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int k = kl; k <= ku; k++) {
                G[IDXG(i,0,k)] = V[IDXV(i,0,k)];
            }
        }
    }
    if (ju == jmax) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int k = kl; k <= ku; k++) {
                G[IDXG(i,jmax,k)] = V[IDXV(i,jmax,k)];
            }
        }
    }

    if (kl == 1) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int j = jl; j <= ju; j++) {
                H[IDXH(i,j,0)] = W[IDXW(i,j,0)];
            }
        }
    }
    if (ku == kmax) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int j = jl; j <= ju; j++) {
                H[IDXH(i,j,kmax)] = W[IDXW(i,j,kmax)];
            }
        }
    }

}

void calculateRsMpi(
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        int il, int iu, int jl, int ju, int kl, int ku,
        Real *F, Real *G, Real *H, Real *RS) {
    #pragma omp parallel for
    for (int i = il; i <= iu; i++) {
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                RS[IDXRS(i, j, k)] = ((F[IDXF(i, j, k)] - F[IDXF(i - 1, j, k)]) / dx +
                        (G[IDXG(i, j, k)] - G[IDXG(i, j - 1, k)]) / dy +
                        (H[IDXH(i, j, k)] - H[IDXH(i, j, k - 1)]) / dz) / dt;
            }
        }
    }
}

void calculateDtMpi(
        Real Re, Real Pr, Real tau,
        Real &dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        int il, int iu, int jl, int ju, int kl, int ku,
        Real *U, Real *V, Real *W,
        bool useTemperature) {
    Real uMaxAbs = Real(0.0), vMaxAbs = Real(0.0), wMaxAbs = Real(0.0);

    // First, compute the maximum absolute velocities in x, y and z direction.
    #pragma omp parallel for reduction(max: uMaxAbs)
    for (int i = il-1; i <= iu; i++) {
        for (int j = jl-1; j <= ju+1; j++) {
            for (int k = kl-1; k <= ku+1; k++) {
                uMaxAbs = std::max(uMaxAbs, std::abs(U[IDXU(i,j,k)]));
            }
        }
    }

    #pragma omp parallel for reduction(max: vMaxAbs)
    for (int i = il-1; i <= iu+1; i++) {
        for (int j = jl-1; j <= ju; j++) {
            for (int k = kl-1; k <= ku+1; k++) {
                vMaxAbs = std::max(vMaxAbs, std::abs(V[IDXV(i,j,k)]));
            }
        }
    }

    #pragma omp parallel for reduction(max: wMaxAbs)
    for (int i = il-1; i <= iu+1; i++) {
        for (int j = jl-1; j <= ju+1; j++) {
            for (int k = kl-1; k <= ku; k++) {
                wMaxAbs = std::max(wMaxAbs, std::abs(W[IDXW(i,j,k)]));
            }
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &uMaxAbs, 1, MPI_REAL_CFD3D, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &vMaxAbs, 1, MPI_REAL_CFD3D, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &wMaxAbs, 1, MPI_REAL_CFD3D, MPI_MAX, MPI_COMM_WORLD);

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

    // Now, use formula (14) from worksheet 1 to compute the time step size.
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

void calculateUvwMpi(
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        int il, int iu, int jl, int ju, int kl, int ku,
        int rankL, int rankR, int rankD, int rankU, int rankB, int rankF, Real *bufSend, Real *bufRecv,
        Real *U, Real *V, Real *W, Real *F, Real *G, Real *H, Real *P, FlagType *Flag) {
    #pragma omp parallel for
    for (int i = il-1; i <= iu; i++) {
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                if(isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i+1,j,k)])){
                    U[IDXU(i, j, k)] = F[IDXF(i, j, k)] - dt / dx * (P[IDXP(i + 1, j, k)] - P[IDXP(i, j, k)]);
                }
            }
        }
    }

    #pragma omp parallel for
    for (int i = il; i <= iu; i++) {
        for (int j = jl-1; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                if(isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i,j+1,k)])){
                    V[IDXV(i, j, k)] = G[IDXG(i, j, k)] - dt / dy * (P[IDXP(i, j + 1, k)] - P[IDXP(i, j, k)]);
                }
            }
        }
    }

    #pragma omp parallel for
    for (int i = il; i <= iu; i++) {
        for (int j = jl; j <= ju; j++) {
            for (int k = kl-1; k <= ku; k++) {
                if(isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i,j,k+1)])){
                    W[IDXW(i, j, k)] = H[IDXH(i, j, k)] - dt / dz * (P[IDXP(i, j, k + 1)] - P[IDXP(i, j, k)]);
                }
            }
        }
    }

    MPI_Status status;
    mpiExchangeUvw(U, V, W, il, iu, jl, ju, kl, ku, rankL, rankR, rankD, rankU, rankB, rankF, bufSend, bufRecv, &status);
}

void calculateTemperatureMpi(
        Real Re, Real Pr, Real alpha,
        Real dt, Real dx, Real dy, Real dz,
        int imax, int jmax, int kmax,
        int il, int iu, int jl, int ju, int kl, int ku,
        int rankL, int rankR, int rankD, int rankU, int rankB, int rankF, Real *bufSend, Real *bufRecv,
        Real *U, Real *V, Real *W, Real *T, Real *T_temp, FlagType *Flag) {
    Real duT_dx, dvT_dy, dwT_dz, d2T_dx2, d2T_dy2, d2T_dz2;

    #pragma omp parallel for private(duT_dx, dvT_dy, dwT_dz, d2T_dx2, d2T_dy2, d2T_dz2)
    for (int i = il; i <= iu; i++) {
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
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
    }

    MPI_Status status;
    mpiExchangeCellData(T, il, iu, jl, ju, kl, ku, rankL, rankR, rankD, rankU, rankB, rankF, bufSend, bufRecv, &status);
}
