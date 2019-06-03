//
// Created by christoph on 03.06.19.
//

#ifndef CFD3D_SORSOLVERCPP_HPP
#define CFD3D_SORSOLVERCPP_HPP

#include "Defines.hpp"

void sorSolverCpp(
        Real omg, Real eps, int itermax,
        Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *P, Real *RS, FlagType *Flag);

#endif //CFD3D_SORSOLVERCPP_HPP
