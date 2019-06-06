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

#ifndef CFD3D_INIT_HPP
#define CFD3D_INIT_HPP

#include "Defines.hpp"

/**
 * Initializes the passed arrays U, V, W, P and T. For fluid cells, UI, VI, WI, PI and TI are used.
 * For obstacle cells (or edges between an obstacle cell and a fluid cell or edges between two obstacle cells), the
 * arrays are initialized with 0 instead.
 * @param UI The initialization value for U.
 * @param VI The initialization value for V.
 * @param WI The initialization value for W.
 * @param PI The initialization value for P.
 * @param TI The initialization value for T.
 * @param imax Number of cells in x direction inside of the domain.
 * @param jmax Number of cells in y direction inside of the domain.
 * @param kmax Number of cells in z direction inside of the domain.
 * @param U The velocities in x direction.
 * @param V The velocities in y direction.
 * @param W The velocities in z direction.
 * @param P The pressure values.
 * @param T The temperature values.
 * @param Flag The flag values (@see Flag.hpp for more information).
 */
void initArrays(Real UI, Real VI, Real WI, Real PI, Real TI, int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *P, Real *T, FlagType *Flag);

#endif //CFD3D_INIT_HPP
