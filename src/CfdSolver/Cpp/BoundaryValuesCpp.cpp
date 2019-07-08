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

#include "BoundaryValuesCpp.hpp"
#include "../Flag.hpp"

void setLeftRightBoundariesCpp(
        Real T_h, Real T_c,
        int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *T,
        FlagType *Flag) {
    #pragma omp parallel for
    for (int j = 1; j <= jmax; j++) {
        for (int k = 1; k <= kmax; k++) {
            // Left wall
            if (isNoSlip(Flag[IDXFLAG(0,j,k)])) {
                U[IDXU(0,j,k)] = 0.0;
                V[IDXV(0,j,k)] = -V[IDXV(1,j,k)];
                W[IDXW(0,j,k)] = -W[IDXW(1,j,k)];
            } else if (isFreeSlip(Flag[IDXFLAG(0,j,k)])) {
                U[IDXU(0,j,k)] = 0.0;
                V[IDXV(0,j,k)] = V[IDXV(1,j,k)];
                W[IDXW(0,j,k)] = W[IDXW(1,j,k)];
            } else if (isOutflow(Flag[IDXFLAG(0,j,k)])) {
                U[IDXU(0,j,k)] = U[IDXU(1,j,k)];
                V[IDXV(0,j,k)] = V[IDXV(1,j,k)];
                W[IDXW(0,j,k)] = W[IDXW(1,j,k)];
            }

            // Right wall
            if (isNoSlip(Flag[IDXFLAG(imax+1,j,k)])) {
                U[IDXU(imax,j,k)] = 0.0;
                V[IDXV(imax+1,j,k)] = -V[IDXV(imax,j,k)];
                W[IDXW(imax+1,j,k)] = -W[IDXW(imax,j,k)];
            } else if (isFreeSlip(Flag[IDXFLAG(imax+1,j,k)])) {
                U[IDXU(imax,j,k)] = 0.0;
                V[IDXV(imax+1,j,k)] = V[IDXV(imax,j,k)];
                W[IDXW(imax+1,j,k)] = W[IDXW(imax,j,k)];
            } else if (isOutflow(Flag[IDXFLAG(imax+1,j,k)])) {
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
}

void setDownUpBoundariesCpp(
        Real T_h, Real T_c,
        int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *T,
        FlagType *Flag) {
    #pragma omp parallel for
    for (int i = 1; i <= imax; i++) {
        for (int k = 1; k <= kmax; k++) {
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
}

void setFrontBackBoundariesCpp(
        Real T_h, Real T_c,
        int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *T,
        FlagType *Flag) {
    #pragma omp parallel for
    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            // Back wall
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

            // Front wall
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

            // Back boundary T
            if (isHot(Flag[IDXFLAG(i,j,0)])) {
                T[IDXT(i,j,0)] = 2 * T_h - T[IDXT(i,j,1)];
            }
            else if (isCold(Flag[IDXFLAG(i,j,0)])) {
                T[IDXT(i,j,0)] = 2 * T_c - T[IDXT(i,j,1)];
            }
            else {
                T[IDXT(i,j,0)] = T[IDXT(i,j,1)];
            }

            // Front boundary T
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
}


void setInternalUBoundariesCpp(
        int imax, int jmax, int kmax,
        Real *U,
        FlagType *Flag) {

    for (int i = 1; i <= imax-1; i++) {
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax; k++) {
                int R_check = 0;
                int L_check = 0;
                int R1_check = 0;
                int L1_check = 0;

                if (!isFluid(Flag[IDXFLAG(i, j, k)])) {
                    if (B_R(Flag[IDXFLAG(i, j, k)])) {
                        U[IDXU(i, j, k)] = Real(0);
                        R_check = 1;
                    }

                    if (B_L(Flag[IDXFLAG(i, j, k)])) {
                        U[IDXU(i - 1, j, k)] = Real(0);
                        L_check = 1;
                    }

                    if (B_U(Flag[IDXFLAG(i, j, k)])) {
                        if (L_check == 0) {
                            U[IDXU(i - 1, j, k)] = -U[IDXU(i - 1, j + 1, k)];
                            L1_check = 1;
                        }
                        if (R_check == 0) {
                            U[IDXU(i, j, k)] = -U[IDXU(i, j + 1, k)];
                            R1_check = 1;
                        }
                    }

                    if (B_D(Flag[IDXFLAG(i, j, k)])) {
                        if (L_check == 0) {
                            U[IDXU(i - 1, j, k)] = -U[IDXU(i - 1, j - 1, k)];
                            L1_check = 1;
                        }
                        if (R_check == 0) {
                            U[IDXU(i, j, k)] = -U[IDXU(i, j - 1, k)];
                            R1_check = 1;
                        }
                    }

                    if (B_B(Flag[IDXFLAG(i, j, k)])) {
                        if (L_check == 0 && L1_check == 0) {
                            U[IDXU(i - 1, j, k)] = -U[IDXU(i - 1, j, k - 1)];
                        }
                        if (R_check == 0 && R1_check == 0) {
                            U[IDXU(i, j, k)] = -U[IDXU(i, j, k - 1)];
                        }
                    }

                    if (B_F(Flag[IDXFLAG(i, j, k)])) {
                        if (L_check == 0 && L1_check == 0) {
                            U[IDXU(i - 1, j, k)] = -U[IDXU(i - 1, j, k + 1)];
                        }
                        if (R_check == 0 && R1_check == 0) {
                            U[IDXU(i, j, k)] = -U[IDXU(i, j, k + 1)];
                        }
                    }
                }
            }
        }
    }     
}

void setInternalVBoundariesCpp(
        int imax, int jmax, int kmax,
        Real *V,
        FlagType *Flag) {
    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax-1; j++) {
            for (int k = 1; k <= kmax; k++) {
                int U_check = 0;
                int D_check = 0;
                int U1_check = 0;
                int D1_check = 0;

                if (!isFluid(Flag[IDXFLAG(i, j, k)])) {
                    if (B_U(Flag[IDXFLAG(i, j, k)])) {
                        V[IDXV(i, j, k)] = Real(0);
                        U_check = 1;
                    }

                    if (B_D(Flag[IDXFLAG(i, j, k)])) {
                        V[IDXV(i, j - 1, k)] = Real(0);
                        D_check = 1;
                    }

                    if (B_R(Flag[IDXFLAG(i, j, k)])) {
                        if (D_check == 0) {
                            V[IDXV(i, j - 1, k)] = -V[IDXV(i + 1, j - 1, k)];
                            D1_check = 0;
                        }
                        if (U_check == 0) {
                            V[IDXV(i, j, k)] = -V[IDXV(i + 1, j, k)];
                            U1_check = 0;
                        }
                    }

                    if (B_L(Flag[IDXFLAG(i, j, k)])) {
                        if (D_check == 0) {
                            V[IDXV(i, j - 1, k)] = -V[IDXV(i - 1, j - 1, k)];
                            D1_check = 0;
                        }
                        if (U_check == 0) {
                            V[IDXV(i, j, k)] = -V[IDXV(i - 1, j, k)];
                            U1_check = 0;
                        }
                    }

                    if (B_B(Flag[IDXFLAG(i, j, k)])) {
                        if (D_check == 0 && D1_check == 0) {
                            V[IDXV(i, j - 1, k)] = -V[IDXV(i, j - 1, k - 1)];
                        }
                        if (U_check == 0 && U1_check == 0) {
                            V[IDXV(i, j, k)] = -V[IDXV(i, j, k - 1)];
                        }
                    }

                    if (B_F(Flag[IDXFLAG(i, j, k)])) {
                        if (D_check == 0 && D1_check == 0) {
                            V[IDXV(i, j - 1, k)] = -V[IDXV(i, j - 1, k + 1)];
                        }
                        if (U_check == 0 && U1_check == 0) {
                            V[IDXV(i, j, k)] = -V[IDXV(i, j, k + 1)];
                        }
                    }
                }
            }
        }
    }   
}

void setInternalWBoundariesCpp(
        int imax, int jmax, int kmax,
        Real *W,
        FlagType *Flag) {
    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax-1; k++) {
                int F_check = 0;
                int B_check = 0;
                int F1_check = 0;
                int B1_check = 0;

                if (!isFluid(Flag[IDXFLAG(i, j, k)])) {
                    if (B_B(Flag[IDXFLAG(i, j, k)])) {
                        W[IDXW(i, j, k - 1)] = Real(0);
                        B_check = 1;
                    }

                    if (B_F(Flag[IDXFLAG(i, j, k)])) {
                        W[IDXW(i, j, k)] = Real(0);
                        F_check = 1;
                    }

                    if (B_R(Flag[IDXFLAG(i, j, k)])) {
                        if (B_check == 0) {
                            W[IDXW(i, j, k - 1)] = -W[IDXW(i + 1, j, k - 1)];
                            B1_check = 1;
                        }
                        if (F_check == 0) {
                            W[IDXW(i, j, k)] = -W[IDXW(i + 1, j, k)];
                            F1_check = 1;
                        }
                    }

                    if (B_L(Flag[IDXFLAG(i, j, k)])) {
                        if (B_check == 0) {
                            W[IDXW(i, j, k - 1)] = -W[IDXW(i - 1, j, k - 1)];
                            B1_check = 1;
                        }
                        if (F_check == 0) {
                            W[IDXW(i, j, k)] = -W[IDXW(i - 1, j, k)];
                            F1_check = 1;
                        }
                    }

                    if (B_U(Flag[IDXFLAG(i, j, k)])) {
                        if (B_check == 0 && B1_check == 0) {
                            W[IDXW(i, j, k - 1)] = -W[IDXW(i, j + 1, k - 1)];
                        }
                        if (F_check == 0 && F1_check == 0) {
                            W[IDXW(i, j, k)] = -W[IDXW(i, j + 1, k)];
                        }
                    }

                    if (B_D(Flag[IDXFLAG(i, j, k)])) {
                        if (B_check == 0 && B1_check == 0) {
                            W[IDXW(i, j, k - 1)] = -W[IDXW(i, j - 1, k - 1)];
                        }
                        if (F_check == 0 && F1_check == 0) {
                            W[IDXW(i, j, k)] = -W[IDXW(i, j - 1, k)];
                        }
                    }
                }
            }     
        }
    }      
}

void setInternalTBoundariesCpp(
        int imax, int jmax, int kmax,
        Real *T,
        FlagType *Flag) {
    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax; k++) {
                int numDirectFlag = 0;
                Real T_temp = Real(0);

                if (!isFluid(Flag[IDXFLAG(i, j, k)])) {
                    if (B_R(Flag[IDXFLAG(i, j, k)])) {
                        T_temp = T[IDXT(i + 1, j, k)];
                        numDirectFlag++;
                    }

                    if (B_L(Flag[IDXFLAG(i, j, k)])) {
                        T_temp = T[IDXT(i - 1, j, k)];
                        numDirectFlag++;
                    }

                    if (B_U(Flag[IDXFLAG(i, j, k)])) {
                        T_temp = T[IDXT(i, j + 1, k)];
                        numDirectFlag++;
                    }

                    if (B_D(Flag[IDXFLAG(i, j, k)])) {
                        T_temp = T[IDXT(i, j - 1, k)];
                        numDirectFlag++;
                    }

                    if (B_B(Flag[IDXFLAG(i, j, k)])) {
                        T_temp = T[IDXT(i, j, k - 1)];
                        numDirectFlag++;
                    }

                    if (B_F(Flag[IDXFLAG(i, j, k)])) {
                        T_temp = T[IDXT(i, j, k + 1)];
                        numDirectFlag++;
                    }

                    T[IDXT(i,j,k)] = T_temp/Real(numDirectFlag);
                }
            }     
        }
    }      
}

void setBoundaryValuesCpp(
        Real T_h, Real T_c,
        int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *T,
        FlagType *Flag) {
    setLeftRightBoundariesCpp(T_h, T_c, imax, jmax, kmax, U, V, W, T, Flag);
    setDownUpBoundariesCpp(T_h, T_c, imax, jmax, kmax, U, V, W, T, Flag);
    setFrontBackBoundariesCpp(T_h, T_c, imax, jmax, kmax, U, V, W, T, Flag);
    setInternalUBoundariesCpp(imax, jmax, kmax, U,Flag);
    setInternalVBoundariesCpp(imax, jmax, kmax, V,Flag);
    setInternalWBoundariesCpp(imax, jmax, kmax, W,Flag);
    setInternalTBoundariesCpp(imax, jmax, kmax, T,Flag);
}


void setBoundaryValuesScenarioSpecificCpp(
        const std::string &scenarioName,
        int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W,
        FlagType *Flag) {
    if (scenarioName == "driven_cavity") {
        #pragma omp parallel for
        for (int i = 1; i <= imax; i++) {
            for (int k = 1; k <= kmax; k++) {
                // Upper wall
                U[IDXU(i,jmax+1,k)] = 2.0 - U[IDXU(i,jmax,k)];
            }
        }
    } else if (scenarioName == "flow_over_step") {
        #pragma omp parallel for
        for (int j = jmax/2+1; j <= jmax; j++) {
            for (int k = 1; k <= kmax; k++) {
                // Left wall
                U[IDXU(0,j,k)] = 1.0;
                V[IDXV(0,j,k)] = 0.0;
                W[IDXW(0,j,k)] = 0.0;
                //U[IDXU(0,j,k)] = 0.0;
                //V[IDXV(0,j,k)] = 2.0 - V[IDXV(1,j,k)];
                //W[IDXW(0,j,k)] = -W[IDXW(1,j,k)];
            }
        }
    } else if (scenarioName == "single_tower") {
        #pragma omp parallel for
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax; k++) {
                // Left wall
                U[IDXU(0,j,k)] = 1.0;
                V[IDXV(0,j,k)] = 0.0;
                W[IDXW(0,j,k)] = 0.0;
            }
        }
    } else if (scenarioName == "terrain_1" || scenarioName == "fuji_san" || scenarioName == "zugspitze") {
        #pragma omp parallel for
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax; k++) {
                // Left wall
                U[IDXU(0,j,k)] = 1.0;
                V[IDXV(0,j,k)] = 0.0;
                W[IDXW(0,j,k)] = 0.0;
            }
        }
    }
}
