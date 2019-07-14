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

#ifndef CFD3D_OPENCLDEFINES_HPP
#define CFD3D_OPENCLDEFINES_HPP

/**
 * The type used for the Flag array. The array stores the type of each cell (i.e., fluid, type of obstacle cell, etc.).
 */
#define FlagType unsigned int


/**
 * The floating point type used for the simulation. float is faster, but double has a higher accuracy.
 */
#ifdef REAL_FLOAT
#define Real float
#endif
#ifdef REAL_DOUBLE
#define Real double
#endif

// Computes the square of a number.
#define SQR(x) ((x)*(x))

/**
 * For accessing 1D-arrays as 3D-arrays.
 */
#define IDXU(i,j,k) ((i)*(jmax+2)*(kmax+2) + (j)*(kmax+2) + (k))
#define IDXV(i,j,k) ((i)*(jmax+1)*(kmax+2) + (j)*(kmax+2) + (k))
#define IDXW(i,j,k) ((i)*(jmax+2)*(kmax+1) + (j)*(kmax+1) + (k))
#define IDXP(i,j,k) ((i)*(jmax+2)*(kmax+2) + (j)*(kmax+2) + (k))
#define IDXT(i,j,k) ((i)*(jmax+2)*(kmax+2) + (j)*(kmax+2) + (k))
#define IDXF(i,j,k) ((i)*(jmax+1)*(kmax+1) + (j)*(kmax+1) + (k))
#define IDXG(i,j,k) ((i)*(jmax+1)*(kmax+1) + (j)*(kmax+1) + (k))
#define IDXH(i,j,k) ((i)*(jmax+1)*(kmax+1) + (j)*(kmax+1) + (k))
#define IDXRS(i,j,k) ((i)*(jmax+1)*(kmax+1) + (j)*(kmax+1) + (k))
#define IDXFLAG(i,j,k) ((i)*(jmax+2)*(kmax+2) + (j)*(kmax+2) + (k))

inline bool isFluid(unsigned int flag) { return (flag >> 0) & 1; }
inline bool isNoSlip(unsigned int flag) { return (flag >> 1) & 1; }
inline bool isFreeSlip(unsigned int flag) { return (flag >> 2) & 1; }
inline bool isOutflow(unsigned int flag) { return (flag >> 3) & 1; }
inline bool isInflow(unsigned int flag) { return (flag >> 4) & 1; }
inline bool B_L(unsigned int flag) { return (flag >> 5) & 1; }
inline bool B_R(unsigned int flag) { return (flag >> 6) & 1; }
inline bool B_D(unsigned int flag) { return (flag >> 7) & 1; }
inline bool B_U(unsigned int flag) { return (flag >> 8) & 1; }
inline bool B_B(unsigned int flag) { return (flag >> 9) & 1; }
inline bool B_F(unsigned int flag) { return (flag >> 10) & 1; }
inline bool isHot(unsigned int flag) { return (flag >> 11) & 1; }
inline bool isCold(unsigned int flag) { return (flag >> 12) & 1; }
inline bool isCoupling(unsigned int flag) { return (flag >> 13) & 1; }

#endif //CFD3D_OPENCLDEFINES_HPP
