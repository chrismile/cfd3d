//
// Created by christoph on 02.06.19.
//

#ifndef CFD3D_DEFINES_HPP
#define CFD3D_DEFINES_HPP

#include <glm/vec3.hpp>

/**
 * The type used for the Flag array. The array stores the type of each cell (i.e., fluid, type of obstacle cell, etc.).
 */
typedef unsigned int FlagType;

/**
 * The floating point type used for the simulation. float is faster, but double has a higher accuracy.
 */
typedef float Real;
//typedef double Real;
typedef glm::vec<3, Real, glm::defaultp> rvec3;

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

#endif //CFD3D_DEFINES_HPP
