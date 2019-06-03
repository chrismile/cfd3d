//
// Created by christoph on 01.06.19.
//

#include "NetCdfWriter.hpp"

void NetCdfWriter::openFile(const std::string &filename,
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz, Real xOrigin, Real yOrigin, Real zOrigin) {
    this->imax = imax;
    this->jmax = jmax;
    this->kmax = kmax;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
    this->xOrigin = xOrigin;
    this->yOrigin = yOrigin;
    this->zOrigin = zOrigin;

    // TODO
}

void NetCdfWriter::writeTimestep(int timeStepNumber, Real time, Real *U, Real *V, Real *W, Real *P, Real *T,
        FlagType *Flag) {
    // TODO
}