//
// Created by christoph on 03.06.19.
//

#ifndef CFD3D_UVWCPP_HPP
#define CFD3D_UVWCPP_HPP

#include "Defines.hpp"

/*
 * Determines the value of F, H and H for computing RS.
 */
void calculateFghCpp(
        Real Re, Real GX, Real GY, Real GZ, Real alpha, Real beta,
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *T, Real *F, Real *G, Real *H, FlagType *Flag);

/*
 * Computes the right hand side of the Pressure Poisson Equation (PPE).
 */
void calculateRsCpp(
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *F, Real *G, Real *H, Real *RS);

/*
 * Determines the maximum time step size. The time step size is restricted according to the CFL theorem.
 */
void calculateDtCpp(
        Real Re, Real Pr, Real tau,
        Real &dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W,
        bool useTemperature);

/*
 * Calculates the new velocity values.
 */
void calculateUvwCpp(
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *F, Real *G, Real *H, Real *P, FlagType *Flag);

/*
 * Calculates the new temperature values.
 */
void calculateTemperatureCpp(
        Real Re, Real Pr, Real alpha,
        Real dt, Real dx, Real dy, Real dz,
        int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *T, Real *T_temp, FlagType *Flag);

#endif //CFD3D_UVWCPP_HPP
