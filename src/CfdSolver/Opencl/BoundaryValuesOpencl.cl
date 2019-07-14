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

kernel void setLeftRightBoundariesOpencl(
        Real T_h, Real T_c,
        int imax, int jmax, int kmax,
        global Real *U, global Real *V, global Real *W, global Real *T,
        global FlagType *Flag) {
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    // Set the boundary values for the pressure on the y-z-planes.
    if (j <= jmax && k <= kmax){
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

kernel void setDownUpBoundariesOpencl(
        Real T_h, Real T_c,
        int imax, int jmax, int kmax,
        global Real *U, global Real *V, global Real *W, global Real *T,
        global FlagType *Flag) {
    int i = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    // Set the boundary values for the pressure on the x-z-planes.
    if (i <= imax && k <= kmax){
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

kernel void setFrontBackBoundariesOpencl(
        Real T_h, Real T_c,
        int imax, int jmax, int kmax,
        global Real *U, global Real *V, global Real *W, global Real *T,
        global FlagType *Flag) {
    int i = get_global_id(1) + 1;
    int j = get_global_id(0) + 1;

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

kernel void setInternalUBoundariesOpencl(
        int imax, int jmax, int kmax,
        global Real *U, global FlagType *Flag) {
    int i = get_global_id(2) + 1;
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    if (i <= imax - 1 && j <= jmax && k <= kmax){
        int R_check = 0;
        int L_check = 0;
        int R1_check = 0;
        int L1_check = 0;

        if (!isFluid(Flag[IDXFLAG(i, j, k)])) {
            if (B_R(Flag[IDXFLAG(i, j, k)])) {
                U[IDXU(i, j, k)] = (Real)(0);
                R_check = 1;
            }

            if (B_L(Flag[IDXFLAG(i, j, k)])) {
                U[IDXU(i - 1, j, k)] = (Real)(0);
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

kernel void setInternalVBoundariesOpencl(
        int imax, int jmax, int kmax,
        global Real *V, global FlagType *Flag) {
    int i = get_global_id(2) + 1;
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    if (i <= imax && j <= jmax - 1 && k <= kmax){
        int U_check = 0;
        int D_check = 0;
        int U1_check = 0;
        int D1_check = 0;

        if (!isFluid(Flag[IDXFLAG(i, j, k)])) {
            if (B_U(Flag[IDXFLAG(i, j, k)])) {
                V[IDXV(i, j, k)] = (Real)(0);
                U_check = 1;
            }

            if (B_D(Flag[IDXFLAG(i, j, k)])) {
                V[IDXV(i, j - 1, k)] = (Real)(0);
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

kernel void setInternalWBoundariesOpencl(
        int imax, int jmax, int kmax,
        global Real *W, global FlagType *Flag) {
    int i = get_global_id(2) + 1;
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    if (i <= imax && j <= jmax && k <= kmax - 1){
        int F_check = 0;
        int B_check = 0;
        int F1_check = 0;
        int B1_check = 0;

        if (!isFluid(Flag[IDXFLAG(i, j, k)])) {
            if (B_B(Flag[IDXFLAG(i, j, k)])) {
                W[IDXW(i, j, k - 1)] = (Real)(0);
                B_check = 1;
            }

            if (B_F(Flag[IDXFLAG(i, j, k)])) {
                W[IDXW(i, j, k)] = (Real)(0);
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

kernel void setInternalTBoundariesOpencl(
        int imax, int jmax, int kmax,
        global Real *T, global FlagType *Flag) {
    int i = get_global_id(2) + 1;
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    if (i <= imax && j <= jmax && k <= kmax){
        int numDirectFlag = 0;
        Real T_temp = (Real)(0);

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

            if (numDirectFlag == 0) {
                T[IDXT(i,j,k)] = 0;
            } else {
                T[IDXT(i,j,k)] = T_temp / (Real)(numDirectFlag);
            }
        }
    }
}



kernel void setDrivenCavityBoundariesOpencl(int imax, int jmax, int kmax,
        global Real *U){
    int i = get_global_id(1);
    int k = get_global_id(0) + 1;

    if (i <= imax && k <= kmax) {
        // Upper wall
        U[IDXU(i,jmax+1,k)] = 2.0 - U[IDXU(i,jmax,k)];
    }
}

kernel void setFlowOverStepBoundariesOpencl(int imax, int jmax, int kmax,
        global Real *U, global Real *V, global Real *W){
    int j = get_global_id(1) + jmax/2+1;
    int k = get_global_id(0) + 1;

    if(j <= jmax && k <= kmax) {
        // Left wall
        U[IDXU(0,j,k)] = 1.0;
        V[IDXV(0,j,k)] = 0.0;
        W[IDXW(0,j,k)] = 0.0;
    }
}

kernel void setSingleTowerBoundariesOpencl(int imax, int jmax, int kmax,
        global Real *U, global Real *V, global Real *W){
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    if(j <= jmax && k <= kmax) {
        // Left wall
        U[IDXU(0,j,k)] = 1.0;
        V[IDXV(0,j,k)] = 0.0;
        W[IDXW(0,j,k)] = 0.0;
    }
}

kernel void setMountainBoundariesOpencl(int imax, int jmax, int kmax,
        global Real *U, global Real *V, global Real *W, global FlagType *Flag){
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    if(j <= jmax && k <= kmax){
        // Left wall
        if (isInflow(Flag[IDXFLAG(0, j, k)])) {
            U[IDXU(0, j, k)] = 1.0;
            V[IDXV(0, j, k)] = 0.0;
            W[IDXW(0, j, k)] = 0.0;
        }
    }
}
