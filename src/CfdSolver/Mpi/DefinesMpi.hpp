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

#ifndef CFD3D_DEFINESMPI_HPP
#define CFD3D_DEFINESMPI_HPP

#ifdef IDXU

#undef IDXU
#undef IDXV
#undef IDXW
#undef IDXP
#undef IDXT
#undef IDXF
#undef IDXG
#undef IDXH
#undef IDXRS
#undef IDXFLAG

#endif

#define IDXU(i,j,k) (((i) - (il-2))*(ju - jl + 3)*(ku - kl + 3) + ((j) - (jl-1))*(ku - kl + 3) + ((k) - (kl-1)))
#define IDXV(i,j,k) (((i) - (il-1))*(ju - jl + 4)*(ku - kl + 3) + ((j) - (jl-2))*(ku - kl + 3) + ((k) - (kl-1)))
#define IDXW(i,j,k) (((i) - (il-1))*(ju - jl + 3)*(ku - kl + 4) + ((j) - (jl-1))*(ku - kl + 4) + ((k) - (kl-2)))
#define IDXP(i,j,k) (((i) - (il-1))*(ju - jl + 3)*(ku - kl + 3) + ((j) - (jl-1))*(ku - kl + 3) + ((k) - (kl-1)))
#define IDXT(i,j,k) (((i) - (il-1))*(ju - jl + 3)*(ku - kl + 3) + ((j) - (jl-1))*(ku - kl + 3) + ((k) - (kl-1)))
#define IDXF(i,j,k) (((i) - (il-2))*(ju - jl + 3)*(ku - kl + 3) + ((j) - (jl-1))*(ku - kl + 3) + ((k) - (kl-1)))
#define IDXG(i,j,k) (((i) - (il-1))*(ju - jl + 4)*(ku - kl + 3) + ((j) - (jl-2))*(ku - kl + 3) + ((k) - (kl-1)))
#define IDXH(i,j,k) (((i) - (il-1))*(ju - jl + 3)*(ku - kl + 4) + ((j) - (jl-1))*(ku - kl + 4) + ((k) - (kl-2)))
#define IDXRS(i,j,k) (((i) - (il))*(ju - jl + 1)*(ku - kl + 1) + ((j) - (jl))*(ku - kl + 1) + ((k) - (kl)))
#define IDXFLAG(i,j,k) (((i) - (il-1))*(ju - jl + 3)*(ku - kl + 3) + ((j) - (jl-1))*(ku - kl + 3) + ((k) - (kl-1)))

#ifdef REAL_FLOAT
#define MPI_REAL_CFD3D MPI_FLOAT
#endif
#ifdef REAL_DOUBLE
#define MPI_REAL_CFD3D MPI_DOUBLE
#endif

#endif //CFD3D_DEFINESMPI_HPP
