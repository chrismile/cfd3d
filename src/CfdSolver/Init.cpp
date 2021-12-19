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

#include "Init.hpp"

void initArrays(
        Real UI, Real VI, Real WI, Real PI, Real TI, int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *P, Real *T, FlagType *Flag) {
    #pragma omp parallel for shared(imax, jmax, kmax, U, UI) default(none)
    for (int i = 0; i <= imax; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            for (int k = 0; k <= kmax+1; k++) {
                U[IDXU(i,j,k)] = UI;
            }
        }
    }

    #pragma omp parallel for shared(imax, jmax, kmax, V, VI) default(none)
    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax; j++) {
            for (int k = 0; k <= kmax+1; k++) {
                V[IDXV(i,j,k)] = VI;
            }
        }
    }

    #pragma omp parallel for shared(imax, jmax, kmax, W, WI) default(none)
    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            for (int k = 0; k <= kmax; k++) {
                W[IDXW(i,j,k)] = WI;
            }
        }
    }

    #pragma omp parallel for shared(imax, jmax, kmax, P, T, PI, TI) default(none)
    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            for (int k = 0; k <= kmax+1; k++) {
                P[IDXP(i,j,k)] = PI;
                T[IDXT(i,j,k)] = TI;
            }
        }
    }
}

#include "CfdSolver/Mpi/DefinesMpi.hpp"

void initArraysMpi(
        Real UI, Real VI, Real WI, Real PI, Real TI, int il, int iu, int jl, int ju, int kl, int ku,
        Real *U, Real *V, Real *W, Real *P, Real *T, FlagType *Flag) {
    for (int i = il - 2; i <= iu + 1; i++) {
        for (int j = jl - 1; j <= ju + 1; j++) {
            for (int k = kl - 1; k <= ku + 1; k++) {
                U[IDXU(i,j,k)] = UI;
            }
        }
    }

    for (int i = il - 1; i <= iu + 1; i++) {
        for (int j = jl - 2; j <= ju + 1; j++) {
            for (int k = kl - 1; k <= ku + 1; k++) {
                V[IDXV(i,j,k)] = VI;
            }
        }
    }

    for (int i = il - 1; i <= iu + 1; i++) {
        for (int j = jl - 1; j <= ju + 1; j++) {
            for (int k = kl - 2; k <= ku + 1; k++) {
                W[IDXW(i,j,k)] = WI;
            }
        }
    }

    for (int i = il - 1; i <= iu + 1; i++) {
        for (int j = jl - 1; j <= ju + 1; j++) {
            for (int k = kl - 1; k <= ku + 1; k++) {
                P[IDXP(i,j,k)] = PI;
                T[IDXT(i,j,k)] = TI;
            }
        }
    }
}
