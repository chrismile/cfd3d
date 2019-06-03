//
// Created by christoph on 01.06.19.
//

#include <cassert>
#include "ParticleTracer.hpp"

Real getUAtIdx(const glm::ivec3 &staggeredGridPosition, int imax, int jmax, int kmax, Real *U) {
    assert(glm::all(glm::greaterThanEqual(staggeredGridPosition, glm::ivec3(0,0,0)))
            && glm::all(glm::lessThan(staggeredGridPosition, glm::ivec3(imax+1,jmax+2,jmax+2))));
    return U[IDXU(staggeredGridPosition.x, staggeredGridPosition.y, staggeredGridPosition.z)];
}
Real getVAtIdx(const glm::ivec3 &staggeredGridPosition, int imax, int jmax, int kmax, Real *V) {
    assert(glm::all(glm::greaterThanEqual(staggeredGridPosition, glm::ivec3(0,0,0)))
           && glm::all(glm::lessThan(staggeredGridPosition, glm::ivec3(imax+2,jmax+1,jmax+2))));
    return V[IDXV(staggeredGridPosition.x, staggeredGridPosition.y, staggeredGridPosition.z)];
}
Real getWAtIdx(const glm::ivec3 &staggeredGridPosition, int imax, int jmax, int kmax, Real *W) {
    assert(glm::all(glm::greaterThanEqual(staggeredGridPosition, glm::ivec3(0,0,0)))
           && glm::all(glm::lessThan(staggeredGridPosition, glm::ivec3(imax+2,jmax+2,jmax+1))));
    return W[IDXW(staggeredGridPosition.x, staggeredGridPosition.y, staggeredGridPosition.z)];
}

Real trilinearInterpolationU(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real *U) {
    rvec3 staggeredGridRealPosition = (particlePosition - gridOrigin) / gridSize + rvec3(0, 0.5, 0.5);
    glm::ivec3 staggeredGridPosition = glm::ivec3(particlePosition);
    rvec3 frac = glm::fract(particlePosition);
    rvec3 invFrac = rvec3(1.0) - frac;
    Real interpolationValue =
            invFrac.x*invFrac.y*invFrac.z*getUAtIdx(staggeredGridPosition + glm::ivec3(0,0,0), imax, kmax, jmax, U)
            + frac.x*invFrac.y*invFrac.z*getUAtIdx(staggeredGridPosition + glm::ivec3(1,0,0), imax, kmax, jmax, U)
            + invFrac.x*frac.y*invFrac.z*getUAtIdx(staggeredGridPosition + glm::ivec3(0,1,0), imax, kmax, jmax, U)
            + frac.x*frac.y*invFrac.z*getUAtIdx(staggeredGridPosition + glm::ivec3(1,1,0), imax, kmax, jmax, U)
            + invFrac.x*invFrac.y*frac.z*getUAtIdx(staggeredGridPosition + glm::ivec3(0,0,1), imax, kmax, jmax, U)
            + frac.x*invFrac.y*frac.z*getUAtIdx(staggeredGridPosition + glm::ivec3(1,0,1), imax, kmax, jmax, U)
            + invFrac.x*frac.y*frac.z*getUAtIdx(staggeredGridPosition + glm::ivec3(0,1,1), imax, kmax, jmax, U)
            + frac.x*frac.y*frac.z*getUAtIdx(staggeredGridPosition + glm::ivec3(1,1,1), imax, kmax, jmax, U);
    return interpolationValue;
}

Real trilinearInterpolationV(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real *V) {
    rvec3 staggeredGridRealPosition = (particlePosition - gridOrigin) / gridSize + rvec3(0.5, 0, 0.5);
    glm::ivec3 staggeredGridPosition = glm::ivec3(particlePosition);
    rvec3 frac = glm::fract(particlePosition);
    rvec3 invFrac = rvec3(1.0) - frac;
    Real interpolationValue =
            invFrac.x*invFrac.y*invFrac.z*getVAtIdx(staggeredGridPosition + glm::ivec3(0,0,0), imax, kmax, jmax, V)
            + frac.x*invFrac.y*invFrac.z*getVAtIdx(staggeredGridPosition + glm::ivec3(1,0,0), imax, kmax, jmax, V)
            + invFrac.x*frac.y*invFrac.z*getVAtIdx(staggeredGridPosition + glm::ivec3(0,1,0), imax, kmax, jmax, V)
            + frac.x*frac.y*invFrac.z*getVAtIdx(staggeredGridPosition + glm::ivec3(1,1,0), imax, kmax, jmax, V)
            + invFrac.x*invFrac.y*frac.z*getVAtIdx(staggeredGridPosition + glm::ivec3(0,0,1), imax, kmax, jmax, V)
            + frac.x*invFrac.y*frac.z*getVAtIdx(staggeredGridPosition + glm::ivec3(1,0,1), imax, kmax, jmax, V)
            + invFrac.x*frac.y*frac.z*getVAtIdx(staggeredGridPosition + glm::ivec3(0,1,1), imax, kmax, jmax, V)
            + frac.x*frac.y*frac.z*getVAtIdx(staggeredGridPosition + glm::ivec3(1,1,1), imax, kmax, jmax, V);
    return interpolationValue;
}

Real trilinearInterpolationW(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real *W) {
    rvec3 staggeredGridRealPosition = (particlePosition - gridOrigin) / gridSize + rvec3(0.5, 0.5, 0);
    glm::ivec3 staggeredGridPosition = glm::ivec3(particlePosition);
    rvec3 frac = glm::fract(particlePosition);
    rvec3 invFrac = rvec3(1.0) - frac;
    Real interpolationValue =
            invFrac.x*invFrac.y*invFrac.z*getWAtIdx(staggeredGridPosition + glm::ivec3(0,0,0), imax, kmax, jmax, W)
            + frac.x*invFrac.y*invFrac.z*getWAtIdx(staggeredGridPosition + glm::ivec3(1,0,0), imax, kmax, jmax, W)
            + invFrac.x*frac.y*invFrac.z*getWAtIdx(staggeredGridPosition + glm::ivec3(0,1,0), imax, kmax, jmax, W)
            + frac.x*frac.y*invFrac.z*getWAtIdx(staggeredGridPosition + glm::ivec3(1,1,0), imax, kmax, jmax, W)
            + invFrac.x*invFrac.y*frac.z*getWAtIdx(staggeredGridPosition + glm::ivec3(0,0,1), imax, kmax, jmax, W)
            + frac.x*invFrac.y*frac.z*getWAtIdx(staggeredGridPosition + glm::ivec3(1,0,1), imax, kmax, jmax, W)
            + invFrac.x*frac.y*frac.z*getWAtIdx(staggeredGridPosition + glm::ivec3(0,1,1), imax, kmax, jmax, W)
            + frac.x*frac.y*frac.z*getWAtIdx(staggeredGridPosition + glm::ivec3(1,1,1), imax, kmax, jmax, W);
    return interpolationValue;
}

rvec3 getVectorVelocityAt(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real *U, Real *V, Real *W) {
    Real u = trilinearInterpolationU(particlePosition, gridOrigin, gridSize, imax, jmax, kmax, U);
    Real v = trilinearInterpolationV(particlePosition, gridOrigin, gridSize, imax, jmax, kmax, V);
    Real w = trilinearInterpolationW(particlePosition, gridOrigin, gridSize, imax, jmax, kmax, W);
    return rvec3(u, v, w);
}

rvec3 integrateParticlePositionEuler(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real *U, Real *V, Real *W, Real dt) {
    rvec3 newPosition = particlePosition + dt*getVectorVelocityAt(
            particlePosition, gridOrigin, gridSize, imax, jmax, kmax, U, V, W);
    return newPosition;
}
