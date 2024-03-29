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

#ifndef CFD3D_UVWMPI_HPP
#define CFD3D_UVWMPI_HPP

#include "Defines.hpp"

/*
 * Determines the value of F, H and H for computing RS.
 */
void calculateFghMpi(
        Real Re, Real GX, Real GY, Real GZ, Real alpha, Real beta,
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        int il, int iu, int jl, int ju, int kl, int ku,
        Real *U, Real *V, Real *W, const Real *T, Real *F, Real *G, Real *H, FlagType *Flag);

/*
 * Computes the right hand side of the Pressure Poisson Equation (PPE).
 */
void calculateRsMpi(
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        int il, int iu, int jl, int ju, int kl, int ku,
        const Real *F, const Real *G, const Real *H, Real *RS);

/*
 * Determines the maximum time step size. The time step size is restricted according to the CFL theorem.
 */
void calculateDtMpi(
        Real Re, Real Pr, Real tau,
        Real &dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        int il, int iu, int jl, int ju, int kl, int ku,
        Real *U, Real *V, Real *W,
        bool useTemperature);

/*
 * Calculates the new velocity values.
 */
void calculateUvwMpi(
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        int il, int iu, int jl, int ju, int kl, int ku,
        int rankL, int rankR, int rankD, int rankU, int rankB, int rankF, Real *bufSend, Real *bufRecv,
        Real *U, Real *V, Real *W, const Real *F, const Real *G, const Real *H, const Real *P, FlagType *Flag);

/*
 * Calculates the new temperature values.
 */
void calculateTemperatureMpi(
        Real Re, Real Pr, Real alpha,
        Real dt, Real dx, Real dy, Real dz,
        int imax, int jmax, int kmax,
        int il, int iu, int jl, int ju, int kl, int ku,
        int rankL, int rankR, int rankD, int rankU, int rankB, int rankF, Real *bufSend, Real *bufRecv,
        Real *U, Real *V, Real *W, Real *T, const Real *T_temp, FlagType *Flag);

#endif //CFD3D_UVWMPI_HPP
