//
// Created by christoph on 03.06.19.
//

#include <iostream>
#include "SorSolverCpp.hpp"

void sorSolverCpp(
        Real omg, Real eps, int itermax,
        Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *P, Real *RS, FlagType *Flag) {
    Real residual = 1e9;
    int it = 0;
    while (it < itermax && residual > eps) {
        // TODO: Implement.

        it++;
    }

    if (residual > eps && it == itermax) {
        std::cerr << "\nSOR solver reached maximum number of iterations without converging." << std::endl;
    }
}
