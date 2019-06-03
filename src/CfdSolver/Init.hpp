//
// Created by christoph on 03.06.19.
//

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
