//
// Created by christoph on 03.06.19.
//

#include "UvwCpp.hpp"

void calculateFghCpp(
        Real Re, Real GX, Real GY, Real GZ, Real alpha, Real beta,
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *T, Real *F, Real *G, Real *H, FlagType *Flag) {
    // TODO
}

void calculateRsCpp(
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *F, Real *G, Real *H, Real *RS) {
    // TODO
}

void calculateDtCpp(
        Real Re, Real Pr, Real tau,
        Real &dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W,
        bool useTemperature) {
    // TODO
}

void calculateUvwCpp(
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *F, Real *G, Real *H, Real *P, FlagType *Flag) {
    // TODO
}

void calculateTemperatureCpp(
        Real Re, Real Pr, Real alpha,
        Real dt, Real dx, Real dy, Real dz,
        int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *T, Real *T_temp, FlagType *Flag) {
    // TODO
}
