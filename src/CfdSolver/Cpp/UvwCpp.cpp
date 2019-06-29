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
#include "UvwCpp.hpp"

void calculateFghCpp(
        Real Re, Real GX, Real GY, Real GZ, Real alpha, Real beta,
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *T, Real *F, Real *G, Real *H, FlagType *Flag) {
    
    Real d2u_dx2,d2u_dy2,d2u_dz2,
         d2v_dx2,d2v_dy2,d2v_dz2,
         d2w_dx2,d2w_dy2,d2w_dz2;
    
    Real du2_dx,duv_dy,duw_dz,
         duv_dx,dv2_dy,dvw_dz,
         duw_dx,dvw_dy,dw2_dz;
    
    Real Dx = 1/dx, Dy = 1/dy, Dz = 1/dz;

    for (int i = 1; i <= imax-1; i++) {
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax; k++) {
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
                
                F[IDXF(i,j,k)] = U[IDXU(i,j,k)] + 
                                dt*(
                                    (1/Re)*(d2u_dx2+d2u_dy2+d2u_dz2)-
                                    du2_dx-duv_dy-duw_dz+
                                    GX-(beta/2)*(T[IDXT(i,j,k)]+T[IDXT(i+1,j,k)])*GX
                                    );

            }
        }
    }

    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax-1; j++) {
            for (int k = 1; k <= kmax; k++) {
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
                
                G[IDXG(i,j,k)] = V[IDXV(i,j,k)] + 
                                dt*(
                                    (1/Re)*(d2v_dx2+d2v_dy2+d2v_dz2)-
                                    duv_dx-dv2_dy-dvw_dz+
                                    GY-(beta/2)*(T[IDXT(i,j,k)]+T[IDXT(i,j+1,k)])*GY
                                    );

            }
        }
    } 

    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax-1; k++) {
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
                
                H[IDXH(i,j,k)] = W[IDXW(i,j,k)] + 
                                dt*(
                                    (1/Re)*(d2w_dx2+d2w_dy2+d2w_dz2)-
                                    duw_dx-dvw_dy-dw2_dz+
                                    GZ-(beta/2)*(T[IDXT(i,j,k)]+T[IDXT(i,j,k+1)])*GZ
                                    );

            }
        }
    }

}

void calculateRsCpp(
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *F, Real *G, Real *H, Real *RS) {

    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax; k++) {
                RS[IDXRS(i, j, k)] = ((F[IDXF(i, j, k)] - F[IDXF(i - 1, j, k)]) / dx +
                        (G[IDXG(i, j, k)] - G[IDXG(i, j - 1, k)]) / dy +
                        (H[IDXH(i, j, k)] - H[IDXH(i, j, k - 1)]) / dz) / dt;
            }
        }
    }
}

void calculateDtCpp(
        Real Re, Real Pr, Real tau,
        Real &dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W,
        bool useTemperature) {
    Real uMaxAbs = Real(0.0), vMaxAbs = Real(0.0), wMaxAbs = Real(0.0);

    // First, compute the maximum absolute velocities in x, y and z direction.
    #pragma omp parallel for reduction(max: uMaxAbs)
    for (int i = 0; i <= imax; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            for (int k = 0; k <= kmax+1; k++) {
                uMaxAbs = std::max(uMaxAbs, std::abs(U[IDXU(i,j,k)]));
            }
        }
    }
    #pragma omp parallel for reduction(max: vMaxAbs)
    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax; j++) {
            for (int k = 0; k <= kmax+1; k++) {
                vMaxAbs = std::max(vMaxAbs, std::abs(V[IDXV(i,j,k)]));
            }
        }
    }

    #pragma omp parallel for reduction(max: wMaxAbs)
    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            for (int k = 0; k <= kmax; k++) {
                wMaxAbs = std::max(wMaxAbs, std::abs(V[IDXW(i,j,k)]));
            }
        }
    }

    if (tau < Real(0.0)) {
        // Constant time step manually specified in configuration file. Check for stability.
        assert(2 / Re * dt < dx * dx * dy * dy / (dx * dx + dy * dy));
        assert(uMaxAbs * dt < dx);
        assert(vMaxAbs * dt < dy);
        if (useTemperature){
            assert(dt < (Re*Pr/2)*(1/((1/(dx*dx))+1/(dy*dy))));
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

void calculateUvwCpp(
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *F, Real *G, Real *H, Real *P, FlagType *Flag) {
 
    for (int i = 1; i <= imax - 1; i++) {
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax; k++) {
                U[IDXU(i, j, k)] = F[IDXF(i, j, k)] - dt / dx * (P[IDXP(i + 1, j, k)] - P[IDXP(i, j, k)]);
            }
        }
    }

    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax - 1; j++) {
            for (int k = 1; k <= kmax; k++) {
                V[IDXV(i, j, k)] = G[IDXG(i, j, k)] - dt / dy * (P[IDXP(i, j + 1, k)] - P[IDXP(i, j, k)]);
            }
        }
    }

    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax - 1; k++) {
                W[IDXW(i, j, k)] = H[IDXH(i, j, k)] - dt / dz * (P[IDXP(i, j, k + 1)] - P[IDXP(i, j, k)]);
            }
        }
    }
}

void calculateTemperatureCpp(
        Real Re, Real Pr, Real alpha,
        Real dt, Real dx, Real dy, Real dz,
        int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *T, Real *T_temp, FlagType *Flag) {

    Real duT_dx, dvT_dy, dwT_dz, d2T_dx2, d2T_dy2, d2T_dz2;
    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax; k++) {
                duT_dx =
                    1 / dx * (
                        U[IDXU(i, j, k)] * ((T_temp[IDXT(i, j, k)] + T_temp[IDXT(i + 1, j, k)]) / 2) -
                        U[IDXU(i - 1, j, k)] * ((T_temp[IDXT(i - 1, j, k)] + T_temp[IDXT(i, j, k)]) / 2) +
                        alpha * (
                            std::abs(U[IDXU(i, j, k)])*((T_temp[IDXT(i, j, k)] - T_temp[IDXT(i + 1, j, k)]) / 2) -
                            std::abs(U[IDXU(i - 1, j, k)])*((T_temp[IDXT(i - 1, j, k)] - T_temp[IDXT(i, j, k)]) / 2)
                            )
                        );

                dvT_dy =
                    1 / dy * (
                        V[IDXV(i, j, k)] * ((T_temp[IDXT(i, j, k)] + T_temp[IDXT(i, j + 1, k)]) / 2) -
                        V[IDXV(i, j - 1, k)] * ((T_temp[IDXT(i, j - 1, k)] + T_temp[IDXT(i, j, k)]) / 2) +
                        alpha * (
                            std::abs(V[IDXV(i, j, k)])*((T_temp[IDXT(i, j, k)] - T_temp[IDXT(i, j + 1, k)]) / 2) -
                            std::abs(V[IDXV(i, j - 1, k)])*((T_temp[IDXT(i, j - 1, k)] - T_temp[IDXT(i, j, k)]) / 2)
                            )
                        );

                dwT_dz =
                    1 / dz * (
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
