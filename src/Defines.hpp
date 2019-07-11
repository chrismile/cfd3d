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

#ifndef CFD3D_DEFINES_HPP
#define CFD3D_DEFINES_HPP

#include <glm/vec3.hpp>

#define REAL_FLOAT

/**
 * The type used for the Flag array. The array stores the type of each cell (i.e., fluid, type of obstacle cell, etc.).
 */
typedef unsigned int FlagType;

/**
 * The solver type for the linear system of equations for the Pressure Poisson Equation (PPE).
 */
enum LinearSystemSolverType {
    LINEAR_SOLVER_JACOBI, LINEAR_SOLVER_SOR, LINEAR_SOLVER_GAUSS_SEIDEL, LINEAR_SOLVER_SOR_PARALLEL,
    LINEAR_SOLVER_GAUSS_SEIDEL_PARALLEL
};


/**
 * The floating point type used for the simulation. float is faster, but double has a higher accuracy.
 */
#ifdef REAL_FLOAT
typedef float Real;
#define stringToReal std::stof
#define nc_put_var_real nc_put_var_float
#define nc_put_vara_real nc_put_vara_float
#define nc_put_var1_real nc_put_var1_float
#define NC_REAL 5
#endif
#ifdef REAL_DOUBLE
typedef double Real;
#define stringToReal std::stod
#define nc_put_var_real nc_put_var_double
#define nc_put_vara_real nc_put_vara_double
#define nc_put_var1_real nc_put_var1_double
#define NC_REAL 6
#endif

typedef glm::vec<3, Real, glm::defaultp> rvec3;

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


inline int iceil(int x, int y) { return (x - 1) / y + 1; }

#define blockSize 4

#endif //CFD3D_DEFINES_HPP
