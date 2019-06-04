//
// Created by christoph on 01.06.19.
//

#include <cassert>
#include "ParticleTracer.hpp"

rvec3 integrateParticlePositionEuler(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real *U, Real *V, Real *W, Real dt) {
    rvec3 newPosition = particlePosition + dt*getVectorVelocityAt(
            particlePosition, gridOrigin, gridSize, imax, jmax, kmax, U, V, W);
    return newPosition;
}
