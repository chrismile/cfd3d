//
// Created by christoph on 01.06.19.
//

#include "ParticleTracer.hpp"

Real trilinearInterpolationU(const rvec3 &particlePosition, int imax, int jmax, int kmax, Real *U) {
    // TODO: Implement.
    return 0;
}

rvec3 getVectorVelocityAt(const rvec3 &particlePosition, int imax, int jmax, int kmax, Real *U, Real *V, Real *W) {
    glm::ivec3 particleGridPosition = glm::ivec3(particlePosition);
    rvec3 fractions = glm::fract(particlePosition);
    rvec3 invFractions = rvec3(1.0) - fractions;

    /*float interpolationU = fractions.x * U[IDXU(particleGridPosition.x, particleGridPosition.y, particleGridPosition.z)]
                           + invFractions.x * U[IDXU(particleGridPosition.x-1, particleGridPosition.y, particleGridPosition.z)];
    float interpolationV = fractions.y * V[IDXV(particleGridPosition.x, particleGridPosition.y, particleGridPosition.z)]
                           + invFractions.y * U[IDXV(particleGridPosition.x, particleGridPosition.y-1, particleGridPosition.z)];
    float interpolationW = fractions.z * W[IDXW(particleGridPosition.x, particleGridPosition.y, particleGridPosition.z)]
                           + invFractions.z * U[IDXW(particleGridPosition.x, particleGridPosition.y, particleGridPosition.z-1)];
    return rvec3(interpolationU, interpolationV, interpolationW);*/

    //Real interpU = ;

    // TODO: Implement.
    return rvec3(0);
}

rvec3 integrateParticlePositionEuler(
        const rvec3 &particlePosition, int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real dt) {
    rvec3 newPosition = particlePosition + dt*getVectorVelocityAt(particlePosition, imax, jmax, kmax, U, V, W);
    return newPosition;
}