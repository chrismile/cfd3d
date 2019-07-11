/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2019, Christoph Neuhauser
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
 * FOR ANY DiuECT, INDiuECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef CFD3D_MPIHELPERS_HPP
#define CFD3D_MPIHELPERS_HPP

//#include "mpi.h"
#include "mpi.h"
#include "Defines.hpp"

void mpiInit(
        int argc, char **argv,
        int iproc, int jproc, int kproc, int imax, int jmax, int kmax,
        int &myrank, int &il, int &iu, int &jl, int &ju, int &kl, int &ku,
        int &rankL, int &rankR, int &rankD, int &rankU, int &rankB, int &rankF,
        int &threadIdxI, int &threadIdxJ, int &threadIdxK, int &nproc);

void mpiExchangeCellData(
        Real *PT, int il, int iu, int jl, int ju, int kl, int ku,
        int rankL, int rankR, int rankD, int rankU, int rankB, int rankF,
        Real *bufSend, Real *bufRecv, MPI_Status *status);

void mpiExchangeUvw(
        Real *U, Real *V, Real *W, int il, int iu, int jl, int ju, int kl, int ku,
        int rankL, int rankR, int rankD, int rankU, int rankB, int rankF,
        Real *bufSend, Real *bufRecv, MPI_Status *status);

void mpiStop();

#endif //CFD3D_MPIHELPERS_HPP
