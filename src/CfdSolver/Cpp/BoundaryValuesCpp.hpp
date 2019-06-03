//
// Created by christoph on 03.06.19.
//

#ifndef CFD3D_BOUNDARYVALUESCPP_HPP
#define CFD3D_BOUNDARYVALUESCPP_HPP

#include <string>
#include "Defines.hpp"

/**
 * Sets the boundary condition values of U, V, W and T using the Flag array.
 */
void setBoundaryValuesCpp(
        Real T_h, Real T_c,
        int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *T,
        FlagType *Flag);

/**
 * Sets special boundary conditions (typically something like inflow) specific to the different scenarios.
 */
void setBoundaryValuesScenarioSpecificCpp(
        const std::string &scenarioName,
        int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W,
        FlagType *Flag);

#endif //CFD3D_BOUNDARYVALUESCPP_HPP
