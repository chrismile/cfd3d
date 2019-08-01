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
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <cmath>
#include "MpiHelpers.hpp"
#include "DefinesMpi.hpp"

void mpiInit(
        int argc, char **argv,
        int iproc, int jproc, int kproc, int imax, int jmax, int kmax,
        int &myrank, int &il, int &iu, int &jl, int &ju, int &kl, int &ku,
        int &rankL, int &rankR, int &rankD, int &rankU, int &rankB, int &rankF,
        int &threadIdxI, int &threadIdxJ, int &threadIdxK, int &nproc)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Compute the number of processes in each dimension.
    assert(nproc == iproc*jproc*kproc);

    // Left neighbor.
    if (myrank % iproc == 0) {
        rankL = MPI_PROC_NULL;
    } else {
        rankL = myrank - 1;
    }

    // Right neighbor.
    if (myrank % iproc == iproc-1) {
        rankR = MPI_PROC_NULL;
    } else {
        rankR = myrank + 1;
    }

    // Bottom ("down") neighbor.
    if ((myrank % (iproc*jproc)) / iproc == 0) {
        rankD = MPI_PROC_NULL;
    } else {
        rankD = myrank - iproc;
    }

    // Top ("up") neighbor.
    if ((myrank % (iproc*jproc)) / iproc == jproc-1) {
        rankU = MPI_PROC_NULL;
    } else {
        rankU = myrank + iproc;
    }

    // Back neighbor.
    if (myrank / (iproc*jproc) == 0) {
        rankB = MPI_PROC_NULL;
    } else {
        rankB = myrank - iproc*jproc;
    }

    // Front neighbor.
    if (myrank / (iproc*jproc) == kproc-1) {
        rankF = MPI_PROC_NULL;
    } else {
        rankF = myrank + iproc*jproc;
    }


    // Index of this rank
    threadIdxI = myrank % iproc;
    threadIdxJ = (myrank % (iproc*jproc)) / iproc;
    threadIdxK = myrank / (iproc*jproc);
}

void mpiDomainDecompositionScheduling(int numElements, int numProcesses, int myrank, int &lower, int &upper) {
    /**
     * Compute fair scheduling for data. Idea:
     * Assuming n is the number of elements, k the number of threads.
     * We want to distribute our n elements to k threads.
     *
     * Assume n mod k != 0. Then we can separate the n elements into two
     * sets (of size lambda1 and lambda2) where each element in lambda1
     * has a = ceil(n/k) elements, and each element in lambda2 has
     * b = floor(n/k) elements.
     *
     * I.e.: n = lambda1 * a + lambda2 * b
     * We choose lambda1 + lambda2 = k, such that each thread gets either
     * a or b elements to process.
     *
     * Therefore: lambda2 = k - lambda1,
     * lambda1 = (n - kb) / (a - b)
     */
    const int k = numProcesses;
    const int n = numElements;
    const int a = n / k;
    const int b = iceil(n, k);
    const int lambda1 = a == b ? 0 : (n - k*b)/(a - b);
    const int lambda2 = k - lambda1;

    if (myrank < lambda1) {
        lower = myrank * a + 1;
        upper = (myrank + 1) * a;
    } else {
        lower = lambda1 * a + (myrank - lambda1) * b + 1;
        upper = lambda1 * a + (myrank - lambda1 + 1) * b;
    }
}

void mpiExchangeCellData(
        Real *PT, int il, int iu, int jl, int ju, int kl, int ku,
        int rankL, int rankR, int rankD, int rankU, int rankB, int rankF,
        Real *bufSend, Real *bufRecv, MPI_Status *status) {
    int chunk;

    // Send to the left, receive from the right.
    if (rankL != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                bufSend[(j - jl) * (ku - kl + 1) + (k - kl)] = PT[IDXP(il,j,k)];
            }
        }
    }
    chunk = (ju - jl + 1) * (ku - kl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankL, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankR, 1, MPI_COMM_WORLD, status);
    if (rankR != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                PT[IDXP(iu + 1,j,k)] = bufRecv[(j - jl) * (ku - kl + 1) + (k - kl)];
            }
        }
    }


    // Send to the right, receive from the left.
    if (rankR != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                bufSend[(j - jl) * (ku - kl + 1) + (k - kl)] = PT[IDXP(iu,j,k)];
            }
        }
    }
    chunk = (ju - jl + 1) * (ku - kl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankR, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankL, 1, MPI_COMM_WORLD, status);
    for (int j = jl; j <= ju && rankL != MPI_PROC_NULL; j++) {
        for (int k = kl; k <= ku; k++) {
            PT[IDXP(il - 1,j,k)] = bufRecv[(j - jl) * (ku - kl + 1) + (k - kl)];
        }
    }
    if (rankL != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                PT[IDXP(il - 1,j,k)] = bufRecv[(j - jl) * (ku - kl + 1) + (k - kl)];
            }
        }
    }


    // Send to the bottom, receive from the top.
    if (rankD != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int k = kl; k <= ku; k++) {
                bufSend[(i - il) * (ku - kl + 1) + (k - kl)] = PT[IDXP(i,jl,k)];
            }
        }
    }
    chunk = (iu - il + 1) * (ku - kl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankD, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankU, 1, MPI_COMM_WORLD, status);
    if (rankU != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int k = kl; k <= ku; k++) {
                PT[IDXP(i,ju + 1,k)] = bufRecv[(i - il) * (ku - kl + 1) + (k - kl)];
            }
        }
    }


    // Send to the top, receive from the bottom.
    if (rankU != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int k = kl; k <= ku; k++) {
                bufSend[(i - il) * (ku - kl + 1) + (k - kl)] = PT[IDXP(i,ju,k)];
            }
        }
    }
    chunk = (iu - il + 1) * (ku - kl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankU, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankD, 1, MPI_COMM_WORLD, status);
    if (rankD != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int k = kl; k <= ku; k++) {
                PT[IDXP(i,jl - 1,k)] = bufRecv[(i - il) * (ku - kl + 1) + (k - kl)];
            }
        }
    }


    // Send to the back, receive from the front.
    if (rankB != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int j = jl; j <= ju; j++) {
                bufSend[(i - il) * (ju - jl + 1) + (j - jl)] = PT[IDXP(i,j,kl)];
            }
        }
    }
    chunk = (iu - il + 1) * (ju - jl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankB, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankF, 1, MPI_COMM_WORLD, status);
    if (rankF != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int j = jl; j <= ju; j++) {
                PT[IDXP(i,j,ku + 1)] = bufRecv[(i - il) * (ju - jl + 1) + (j - jl)];
            }
        }
    }


    // Send to the front, receive from the back.
    if (rankF != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int j = jl; j <= ju; j++) {
                bufSend[(i - il) * (ju - jl + 1) + (j - jl)] = PT[IDXP(i,j,ku)];
            }
        }
    }
    chunk = (iu - il + 1) * (ju - jl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankF, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankB, 1, MPI_COMM_WORLD, status);
    if (rankB != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int j = jl; j <= ju; j++) {
                PT[IDXP(i,j,kl - 1)] = bufRecv[(i - il) * (ju - jl + 1) + (j - jl)];
            }
        }
    }
}

void mpiExchangeUvw(
        Real *U, Real *V, Real *W, int il, int iu, int jl, int ju, int kl, int ku,
        int rankL, int rankR, int rankD, int rankU, int rankB, int rankF,
        Real *bufSend, Real *bufRecv, MPI_Status *status) {
    int chunk;

    // Send to the left, receive from the right (in the order U, V, W).
    if (rankL != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                bufSend[(j - jl) * (ku - kl + 1) + (k - kl)] = U[IDXU(il,j,k)];
            }
        }
    }
    chunk = (ju - jl + 1) * (ku - kl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankL, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankR, 1, MPI_COMM_WORLD, status);
    if (rankR != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                U[IDXU(iu + 1,j,k)] = bufRecv[(j - jl) * (ku - kl + 1) + (k - kl)];
            }
        }
    }

    if (rankL != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl - 1; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                bufSend[(j - jl + 1) * (ku - kl + 1) + (k - kl)] = V[IDXV(il,j,k)];
            }
        }
    }
    chunk = (ju - jl + 2) * (ku - kl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankL, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankR, 1, MPI_COMM_WORLD, status);
    if (rankR != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl - 1; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                V[IDXV(iu + 1,j,k)] = bufRecv[(j - jl + 1) * (ku - kl + 1) + (k - kl)];
            }
        }
    }

    if (rankL != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl; j <= ju; j++) {
            for (int k = kl - 1; k <= ku; k++) {
                bufSend[(j - jl) * (ku - kl + 2) + (k - kl + 1)] = W[IDXW(il,j,k)];
            }
        }
    }
    chunk = (ju - jl + 1) * (ku - kl + 2);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankL, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankR, 1, MPI_COMM_WORLD, status);
    if (rankR != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl; j <= ju; j++) {
            for (int k = kl - 1; k <= ku; k++) {
                W[IDXW(iu + 1,j,k)] = bufRecv[(j - jl) * (ku - kl + 2) + (k - kl + 1)];
            }
        }
    }



    // Send to the right, receive from the left (in the order U, V, W).
    if (rankR != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                bufSend[(j - jl) * (ku - kl + 1) + (k - kl)] = U[IDXU(iu - 1,j,k)];
            }
        }
    }
    chunk = (ju - jl + 1) * (ku - kl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankR, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankL, 1, MPI_COMM_WORLD, status);
    if (rankL != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                U[IDXU(il - 2,j,k)] = bufRecv[(j - jl) * (ku - kl + 1) + (k - kl)];
            }
        }
    }

    if (rankR != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl - 1; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                bufSend[(j - jl + 1) * (ku - kl + 1) + (k - kl)] = V[IDXV(iu,j,k)];
            }
        }
    }
    chunk = (ju - jl + 2) * (ku - kl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankR, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankL, 1, MPI_COMM_WORLD, status);
    if (rankL != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl - 1; j <= ju; j++) {
            for (int k = kl; k <= ku; k++) {
                V[IDXV(il - 1,j,k)] = bufRecv[(j - jl + 1) * (ku - kl + 1) + (k - kl)];
            }
        }
    }

    if (rankR != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl; j <= ju; j++) {
            for (int k = kl - 1; k <= ku; k++) {
                bufSend[(j - jl) * (ku - kl + 2) + (k - kl + 1)] = W[IDXW(iu,j,k)];
            }
        }
    }
    chunk = (ju - jl + 1) * (ku - kl + 2);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankR, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankL, 1, MPI_COMM_WORLD, status);
    if (rankL != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int j = jl; j <= ju; j++) {
            for (int k = kl - 1; k <= ku; k++) {
                W[IDXW(il - 1,j,k)] = bufRecv[(j - jl) * (ku - kl + 2) + (k - kl + 1)];
            }
        }
    }



    // Send to the bottom, receive from the top (in the order U, V, W).
    if (rankD != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il - 1; i <= iu; i++) {
            for (int k = kl; k <= ku; k++) {
                bufSend[(i - il + 1) * (ku - kl + 1) + (k - kl)] = U[IDXU(i,jl,k)];
            }
        }
    }
    chunk = (iu - il + 2) * (ku - kl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankD, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankU, 1, MPI_COMM_WORLD, status);
    for (int i = il - 1; i <= iu && rankU != MPI_PROC_NULL; i++) {
        for (int k = kl; k <= ku; k++) {
            U[IDXU(i,ju + 1,k)] = bufRecv[(i - il + 1) * (ku - kl + 1) + (k - kl)];
        }
    }
    if (rankU != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il - 1; i <= iu; i++) {
            for (int k = kl; k <= ku; k++) {
                U[IDXU(i,ju + 1,k)] = bufRecv[(i - il + 1) * (ku - kl + 1) + (k - kl)];
            }
        }
    }

    if (rankD != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int k = kl; k <= ku; k++) {
                bufSend[(i - il) * (ku - kl + 1) + (k - kl)] = V[IDXV(i,jl,k)];
            }
        }
    }
    chunk = (iu - il + 1) * (ku - kl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankD, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankU, 1, MPI_COMM_WORLD, status);
    if (rankU != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int k = kl; k <= ku; k++) {
                V[IDXV(i,ju + 1,k)] = bufRecv[(i - il) * (ku - kl + 1) + (k - kl)];
            }
        }
    }

    if (rankD != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int k = kl - 1; k <= ku; k++) {
                bufSend[(i - il) * (ku - kl + 2) + (k - kl + 1)] = W[IDXW(i,jl,k)];
            }
        }
    }
    chunk = (iu - il + 1) * (ku - kl + 2);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankD, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankU, 1, MPI_COMM_WORLD, status);
    if (rankU != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int k = kl - 1; k <= ku; k++) {
                W[IDXW(i,ju + 1,k)] = bufRecv[(i - il) * (ku - kl + 2) + (k - kl + 1)];
            }
        }
    }



    // Send to the top, receive from the bottom (in the order U, V, W).
    if (rankU != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il - 1; i <= iu; i++) {
            for (int k = kl; k <= ku; k++) {
                bufSend[(i - il + 1) * (ku - kl + 1) + (k - kl)] = U[IDXU(i,ju,k)];
            }
        }
    }
    chunk = (iu - il + 2) * (ku - kl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankU, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankD, 1, MPI_COMM_WORLD, status);
    if (rankD != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il - 1; i <= iu; i++) {
            for (int k = kl; k <= ku; k++) {
                U[IDXU(i,jl - 1,k)] = bufRecv[(i - il + 1) * (ku - kl + 1) + (k - kl)];
            }
        }
    }

    if (rankU != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int k = kl; k <= ku; k++) {
                bufSend[(i - il) * (ku - kl + 1) + (k - kl)] = V[IDXV(i,ju - 1,k)];
            }
        }
    }
    chunk = (iu - il + 1) * (ku - kl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankU, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankD, 1, MPI_COMM_WORLD, status);
    if (rankD != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int k = kl; k <= ku; k++) {
                V[IDXV(i,jl - 2,k)] = bufRecv[(i - il) * (ku - kl + 1) + (k - kl)];
            }
        }
    }

    if (rankU != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int k = kl - 1; k <= ku; k++) {
                bufSend[(i - il) * (ku - kl + 2) + (k - kl + 1)] = W[IDXW(i,ju,k)];
            }
        }
    }
    chunk = (iu - il + 1) * (ku - kl + 2);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankU, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankD, 1, MPI_COMM_WORLD, status);
    if (rankD != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int k = kl - 1; k <= ku; k++) {
                W[IDXW(i,jl - 1,k)] = bufRecv[(i - il) * (ku - kl + 2) + (k - kl + 1)];
            }
        }
    }



    // Send to the back, receive from the front (in the order U, V, W).
    if (rankB != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il - 1; i <= iu; i++) {
            for (int j = jl; j <= ju; j++) {
                bufSend[(i - il + 1) * (ju - jl + 1) + (j - jl)] = U[IDXU(i,j,kl)];
            }
        }
    }
    chunk = (iu - il + 2) * (ju - jl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankB, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankF, 1, MPI_COMM_WORLD, status);
    if (rankF != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il - 1; i <= iu; i++) {
            for (int j = jl; j <= ju; j++) {
                U[IDXU(i,j,ku + 1)] = bufRecv[(i - il + 1) * (ju - jl + 1) + (j - jl)];
            }
        }
    }

    if (rankB != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int j = jl - 1; j <= ju; j++) {
                bufSend[(i - il) * (ju - jl + 2) + (j - jl + 1)] = V[IDXV(i,j,kl)];
            }
        }
    }
    chunk = (iu - il + 1) * (ju - jl + 2);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankB, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankF, 1, MPI_COMM_WORLD, status);
    if (rankF != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int j = jl - 1; j <= ju; j++) {
                V[IDXV(i,j,ku + 1)] = bufRecv[(i - il) * (ju - jl + 2) + (j - jl + 1)];
            }
        }
    }

    if (rankB != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int j = jl; j <= ju; j++) {
                bufSend[(i - il) * (ju - jl + 1) + (j - jl)] = W[IDXW(i,j,kl)];
            }
        }
    }
    chunk = (iu - il + 1) * (ju - jl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankB, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankF, 1, MPI_COMM_WORLD, status);
    if (rankF != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int j = jl; j <= ju; j++) {
                W[IDXW(i,j,ku + 1)] = bufRecv[(i - il) * (ju - jl + 1) + (j - jl)];
            }
        }
    }



    // Send to the front, receive from the back (in the order U, V, W).
    if (rankF != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il - 1; i <= iu; i++) {
            for (int j = jl; j <= ju; j++) {
                bufSend[(i - il + 1) * (ju - jl + 1) + (j - jl)] = U[IDXU(i,j,ku)];
            }
        }
    }
    chunk = (iu - il + 2) * (ju - jl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankF, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankB, 1, MPI_COMM_WORLD, status);
    if (rankB != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il - 1; i <= iu; i++) {
            for (int j = jl; j <= ju; j++) {
                U[IDXU(i,j,kl - 1)] = bufRecv[(i - il + 1) * (ju - jl + 1) + (j - jl)];
            }
        }
    }

    if (rankF != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int j = jl - 1; j <= ju; j++) {
                bufSend[(i - il) * (ju - jl + 2) + (j - jl + 1)] = V[IDXV(i,j,ku)];
            }
        }
    }
    chunk = (iu - il + 1) * (ju - jl + 2);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankF, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankB, 1, MPI_COMM_WORLD, status);
    if (rankB != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int j = jl - 1; j <= ju; j++) {
                V[IDXV(i,j,kl - 1)] = bufRecv[(i - il) * (ju - jl + 2) + (j - jl + 1)];
            }
        }
    }

    if (rankF != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int j = jl; j <= ju; j++) {
                bufSend[(i - il) * (ju - jl + 1) + (j - jl)] = W[IDXW(i,j,ku - 1)];
            }
        }
    }
    chunk = (iu - il + 1) * (ju - jl + 1);
    MPI_Sendrecv(bufSend, chunk, MPI_REAL_CFD3D, rankF, 1, bufRecv, chunk, MPI_REAL_CFD3D, rankB, 1, MPI_COMM_WORLD, status);
    if (rankB != MPI_PROC_NULL) {
        #pragma omp parallel for
        for (int i = il; i <= iu; i++) {
            for (int j = jl; j <= ju; j++) {
                W[IDXW(i,j,kl - 2)] = bufRecv[(i - il) * (ju - jl + 1) + (j - jl)];
            }
        }
    }
}

void mpiStop() {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
