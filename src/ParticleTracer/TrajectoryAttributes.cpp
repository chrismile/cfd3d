//
// Created by christoph on 05.06.19.
//

#include "TrajectoryAttributes.hpp"

void pushTrajectoryAttributes(
        Trajectory &currentTrajectory, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz, Real *U, Real *V, Real *W, Real *P, Real *T) {
    // Necessary to initialize the data first?
    if (currentTrajectory.attributes.size() == 0) {
        currentTrajectory.attributes.resize(2);
    }

    // First attribute: Curl length (also called vorticity).
    rvec3 particlePosition = currentTrajectory.positions.back();
    currentTrajectory.attributes.at(0).push_back(glm::length(getCurlAt(
            particlePosition, gridOrigin, gridSize, imax, jmax, kmax, dx, dy, dz, U, V, W)));
    // Second attribute: Velocity magnitude.
    currentTrajectory.attributes.at(1).push_back(glm::length(getVectorVelocityAt(
            particlePosition, gridOrigin, gridSize, imax, jmax, kmax, U, V, W)));

    // TODO: Possible future attributes: Pressure at the point, temperature at the point, ...
}

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
    rvec3 frac = glm::fract(staggeredGridRealPosition);
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
    rvec3 frac = glm::fract(staggeredGridRealPosition);
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
    rvec3 frac = glm::fract(staggeredGridRealPosition);
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

Real getdUdyAtIdx(const glm::ivec3 &staggeredGridPosition, int imax, int jmax, int kmax, Real *U, int dy) {
    assert(glm::all(glm::greaterThanEqual(staggeredGridPosition, glm::ivec3(0,0,0)))
           && glm::all(glm::lessThan(staggeredGridPosition, glm::ivec3(imax+1,jmax+1,jmax+2))));
    return (U[IDXU(staggeredGridPosition.x, staggeredGridPosition.y, staggeredGridPosition.z)]
            - U[IDXU(staggeredGridPosition.x, staggeredGridPosition.y-1, staggeredGridPosition.z)])/dy;
}
Real getdUdzAtIdx(const glm::ivec3 &staggeredGridPosition, int imax, int jmax, int kmax, Real *U, int dz) {
    assert(glm::all(glm::greaterThanEqual(staggeredGridPosition, glm::ivec3(0,0,0)))
           && glm::all(glm::lessThan(staggeredGridPosition, glm::ivec3(imax+1,jmax+2,jmax+1))));
    return (U[IDXU(staggeredGridPosition.x, staggeredGridPosition.y, staggeredGridPosition.z)]
            - U[IDXU(staggeredGridPosition.x, staggeredGridPosition.y, staggeredGridPosition.z-1)])/dz;
}
Real getdVdxAtIdx(const glm::ivec3 &staggeredGridPosition, int imax, int jmax, int kmax, Real *V, int dx) {
    assert(glm::all(glm::greaterThanEqual(staggeredGridPosition, glm::ivec3(0,0,0)))
           && glm::all(glm::lessThan(staggeredGridPosition, glm::ivec3(imax+1,jmax+1,jmax+2))));
    return (V[IDXU(staggeredGridPosition.x, staggeredGridPosition.y, staggeredGridPosition.z)]
            - V[IDXU(staggeredGridPosition.x-1, staggeredGridPosition.y, staggeredGridPosition.z)])/dx;
}
Real getdVdzAtIdx(const glm::ivec3 &staggeredGridPosition, int imax, int jmax, int kmax, Real *V, int dz) {
    assert(glm::all(glm::greaterThanEqual(staggeredGridPosition, glm::ivec3(0,0,0)))
           && glm::all(glm::lessThan(staggeredGridPosition, glm::ivec3(imax+2,jmax+1,jmax+1))));
    return (V[IDXV(staggeredGridPosition.x, staggeredGridPosition.y, staggeredGridPosition.z)]
            - V[IDXV(staggeredGridPosition.x, staggeredGridPosition.y, staggeredGridPosition.z-1)])/dz;
}
Real getdWdxAtIdx(const glm::ivec3 &staggeredGridPosition, int imax, int jmax, int kmax, Real *W, int dx) {
    assert(glm::all(glm::greaterThanEqual(staggeredGridPosition, glm::ivec3(0,0,0)))
           && glm::all(glm::lessThan(staggeredGridPosition, glm::ivec3(imax+1,jmax+2,jmax+1))));
    return (W[IDXW(staggeredGridPosition.x, staggeredGridPosition.y, staggeredGridPosition.z)]
            - W[IDXW(staggeredGridPosition.x-1, staggeredGridPosition.y, staggeredGridPosition.z)])/dx;
}
Real getdWdyAtIdx(const glm::ivec3 &staggeredGridPosition, int imax, int jmax, int kmax, Real *W, int dy) {
    assert(glm::all(glm::greaterThanEqual(staggeredGridPosition, glm::ivec3(0,0,0)))
           && glm::all(glm::lessThan(staggeredGridPosition, glm::ivec3(imax+2,jmax+1,jmax+1))));
    return (W[IDXW(staggeredGridPosition.x, staggeredGridPosition.y, staggeredGridPosition.z)]
            - W[IDXW(staggeredGridPosition.x, staggeredGridPosition.y-1, staggeredGridPosition.z)])/dy;
}

Real trilinearInterpolation_dUdy(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real dy, Real *U) {
    rvec3 staggeredGridRealPosition = (particlePosition - gridOrigin) / gridSize + rvec3(0, 0, 0.5);
    glm::ivec3 staggeredGridPosition = glm::ivec3(particlePosition);
    rvec3 frac = glm::fract(staggeredGridRealPosition);
    rvec3 invFrac = rvec3(1.0) - frac;
    Real interpolationValue =
            invFrac.x*invFrac.y*invFrac.z*getdUdyAtIdx(staggeredGridPosition + glm::ivec3(0,0,0), imax, kmax, jmax, U, dy)
            + frac.x*invFrac.y*invFrac.z*getdUdyAtIdx(staggeredGridPosition + glm::ivec3(1,0,0), imax, kmax, jmax, U, dy)
            + invFrac.x*frac.y*invFrac.z*getdUdyAtIdx(staggeredGridPosition + glm::ivec3(0,1,0), imax, kmax, jmax, U, dy)
            + frac.x*frac.y*invFrac.z*getdUdyAtIdx(staggeredGridPosition + glm::ivec3(1,1,0), imax, kmax, jmax, U, dy)
            + invFrac.x*invFrac.y*frac.z*getdUdyAtIdx(staggeredGridPosition + glm::ivec3(0,0,1), imax, kmax, jmax, U, dy)
            + frac.x*invFrac.y*frac.z*getdUdyAtIdx(staggeredGridPosition + glm::ivec3(1,0,1), imax, kmax, jmax, U, dy)
            + invFrac.x*frac.y*frac.z*getdUdyAtIdx(staggeredGridPosition + glm::ivec3(0,1,1), imax, kmax, jmax, U, dy)
            + frac.x*frac.y*frac.z*getdUdyAtIdx(staggeredGridPosition + glm::ivec3(1,1,1), imax, kmax, jmax, U, dy);
    return interpolationValue;
}
Real trilinearInterpolation_dUdz(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real dz, Real *U) {
    rvec3 staggeredGridRealPosition = (particlePosition - gridOrigin) / gridSize + rvec3(0, 0.5, 0);
    glm::ivec3 staggeredGridPosition = glm::ivec3(particlePosition);
    rvec3 frac = glm::fract(staggeredGridRealPosition);
    rvec3 invFrac = rvec3(1.0) - frac;
    Real interpolationValue =
            invFrac.x*invFrac.y*invFrac.z*getdUdzAtIdx(staggeredGridPosition + glm::ivec3(0,0,0), imax, kmax, jmax, U, dz)
            + frac.x*invFrac.y*invFrac.z*getdUdzAtIdx(staggeredGridPosition + glm::ivec3(1,0,0), imax, kmax, jmax, U, dz)
            + invFrac.x*frac.y*invFrac.z*getdUdzAtIdx(staggeredGridPosition + glm::ivec3(0,1,0), imax, kmax, jmax, U, dz)
            + frac.x*frac.y*invFrac.z*getdUdzAtIdx(staggeredGridPosition + glm::ivec3(1,1,0), imax, kmax, jmax, U, dz)
            + invFrac.x*invFrac.y*frac.z*getdUdzAtIdx(staggeredGridPosition + glm::ivec3(0,0,1), imax, kmax, jmax, U, dz)
            + frac.x*invFrac.y*frac.z*getdUdzAtIdx(staggeredGridPosition + glm::ivec3(1,0,1), imax, kmax, jmax, U, dz)
            + invFrac.x*frac.y*frac.z*getdUdzAtIdx(staggeredGridPosition + glm::ivec3(0,1,1), imax, kmax, jmax, U, dz)
            + frac.x*frac.y*frac.z*getdUdzAtIdx(staggeredGridPosition + glm::ivec3(1,1,1), imax, kmax, jmax, U, dz);
    return interpolationValue;
}
Real trilinearInterpolation_dVdx(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real dx, Real *U) {
    rvec3 staggeredGridRealPosition = (particlePosition - gridOrigin) / gridSize + rvec3(0, 0, 0.5);
    glm::ivec3 staggeredGridPosition = glm::ivec3(particlePosition);
    rvec3 frac = glm::fract(staggeredGridRealPosition);
    rvec3 invFrac = rvec3(1.0) - frac;
    Real interpolationValue =
            invFrac.x*invFrac.y*invFrac.z*getdVdxAtIdx(staggeredGridPosition + glm::ivec3(0,0,0), imax, kmax, jmax, U, dx)
            + frac.x*invFrac.y*invFrac.z*getdVdxAtIdx(staggeredGridPosition + glm::ivec3(1,0,0), imax, kmax, jmax, U, dx)
            + invFrac.x*frac.y*invFrac.z*getdVdxAtIdx(staggeredGridPosition + glm::ivec3(0,1,0), imax, kmax, jmax, U, dx)
            + frac.x*frac.y*invFrac.z*getdVdxAtIdx(staggeredGridPosition + glm::ivec3(1,1,0), imax, kmax, jmax, U, dx)
            + invFrac.x*invFrac.y*frac.z*getdVdxAtIdx(staggeredGridPosition + glm::ivec3(0,0,1), imax, kmax, jmax, U, dx)
            + frac.x*invFrac.y*frac.z*getdVdxAtIdx(staggeredGridPosition + glm::ivec3(1,0,1), imax, kmax, jmax, U, dx)
            + invFrac.x*frac.y*frac.z*getdVdxAtIdx(staggeredGridPosition + glm::ivec3(0,1,1), imax, kmax, jmax, U, dx)
            + frac.x*frac.y*frac.z*getdVdxAtIdx(staggeredGridPosition + glm::ivec3(1,1,1), imax, kmax, jmax, U, dx);
    return interpolationValue;
}
Real trilinearInterpolation_dVdz(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real dz, Real *U) {
    rvec3 staggeredGridRealPosition = (particlePosition - gridOrigin) / gridSize + rvec3(0.5, 0, 0);
    glm::ivec3 staggeredGridPosition = glm::ivec3(particlePosition);
    rvec3 frac = glm::fract(staggeredGridRealPosition);
    rvec3 invFrac = rvec3(1.0) - frac;
    Real interpolationValue =
            invFrac.x*invFrac.y*invFrac.z*getdVdzAtIdx(staggeredGridPosition + glm::ivec3(0,0,0), imax, kmax, jmax, U, dz)
            + frac.x*invFrac.y*invFrac.z*getdVdzAtIdx(staggeredGridPosition + glm::ivec3(1,0,0), imax, kmax, jmax, U, dz)
            + invFrac.x*frac.y*invFrac.z*getdVdzAtIdx(staggeredGridPosition + glm::ivec3(0,1,0), imax, kmax, jmax, U, dz)
            + frac.x*frac.y*invFrac.z*getdVdzAtIdx(staggeredGridPosition + glm::ivec3(1,1,0), imax, kmax, jmax, U, dz)
            + invFrac.x*invFrac.y*frac.z*getdVdzAtIdx(staggeredGridPosition + glm::ivec3(0,0,1), imax, kmax, jmax, U, dz)
            + frac.x*invFrac.y*frac.z*getdVdzAtIdx(staggeredGridPosition + glm::ivec3(1,0,1), imax, kmax, jmax, U, dz)
            + invFrac.x*frac.y*frac.z*getdVdzAtIdx(staggeredGridPosition + glm::ivec3(0,1,1), imax, kmax, jmax, U, dz)
            + frac.x*frac.y*frac.z*getdVdzAtIdx(staggeredGridPosition + glm::ivec3(1,1,1), imax, kmax, jmax, U, dz);
    return interpolationValue;
}
Real trilinearInterpolation_dWdx(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real dx, Real *U) {
    rvec3 staggeredGridRealPosition = (particlePosition - gridOrigin) / gridSize + rvec3(0, 0.5, 0);
    glm::ivec3 staggeredGridPosition = glm::ivec3(particlePosition);
    rvec3 frac = glm::fract(staggeredGridRealPosition);
    rvec3 invFrac = rvec3(1.0) - frac;
    Real interpolationValue =
            invFrac.x*invFrac.y*invFrac.z*getdWdxAtIdx(staggeredGridPosition + glm::ivec3(0,0,0), imax, kmax, jmax, U, dx)
            + frac.x*invFrac.y*invFrac.z*getdWdxAtIdx(staggeredGridPosition + glm::ivec3(1,0,0), imax, kmax, jmax, U, dx)
            + invFrac.x*frac.y*invFrac.z*getdWdxAtIdx(staggeredGridPosition + glm::ivec3(0,1,0), imax, kmax, jmax, U, dx)
            + frac.x*frac.y*invFrac.z*getdWdxAtIdx(staggeredGridPosition + glm::ivec3(1,1,0), imax, kmax, jmax, U, dx)
            + invFrac.x*invFrac.y*frac.z*getdWdxAtIdx(staggeredGridPosition + glm::ivec3(0,0,1), imax, kmax, jmax, U, dx)
            + frac.x*invFrac.y*frac.z*getdWdxAtIdx(staggeredGridPosition + glm::ivec3(1,0,1), imax, kmax, jmax, U, dx)
            + invFrac.x*frac.y*frac.z*getdWdxAtIdx(staggeredGridPosition + glm::ivec3(0,1,1), imax, kmax, jmax, U, dx)
            + frac.x*frac.y*frac.z*getdWdxAtIdx(staggeredGridPosition + glm::ivec3(1,1,1), imax, kmax, jmax, U, dx);
    return interpolationValue;
}
Real trilinearInterpolation_dWdy(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real dy, Real *U) {
    rvec3 staggeredGridRealPosition = (particlePosition - gridOrigin) / gridSize + rvec3(0, 0.5, 0);
    glm::ivec3 staggeredGridPosition = glm::ivec3(particlePosition);
    rvec3 frac = glm::fract(staggeredGridRealPosition);
    rvec3 invFrac = rvec3(1.0) - frac;
    Real interpolationValue =
            invFrac.x*invFrac.y*invFrac.z*getdWdyAtIdx(staggeredGridPosition + glm::ivec3(0,0,0), imax, kmax, jmax, U, dy)
            + frac.x*invFrac.y*invFrac.z*getdWdyAtIdx(staggeredGridPosition + glm::ivec3(1,0,0), imax, kmax, jmax, U, dy)
            + invFrac.x*frac.y*invFrac.z*getdWdyAtIdx(staggeredGridPosition + glm::ivec3(0,1,0), imax, kmax, jmax, U, dy)
            + frac.x*frac.y*invFrac.z*getdWdyAtIdx(staggeredGridPosition + glm::ivec3(1,1,0), imax, kmax, jmax, U, dy)
            + invFrac.x*invFrac.y*frac.z*getdWdyAtIdx(staggeredGridPosition + glm::ivec3(0,0,1), imax, kmax, jmax, U, dy)
            + frac.x*invFrac.y*frac.z*getdWdyAtIdx(staggeredGridPosition + glm::ivec3(1,0,1), imax, kmax, jmax, U, dy)
            + invFrac.x*frac.y*frac.z*getdWdyAtIdx(staggeredGridPosition + glm::ivec3(0,1,1), imax, kmax, jmax, U, dy)
            + frac.x*frac.y*frac.z*getdWdyAtIdx(staggeredGridPosition + glm::ivec3(1,1,1), imax, kmax, jmax, U, dy);
    return interpolationValue;
}

rvec3 getCurlAt(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz, Real *U, Real *V, Real *W) {
    Real dUdy = trilinearInterpolation_dUdy(particlePosition, gridOrigin, gridSize, imax, jmax, kmax, dy, U);
    Real dUdz = trilinearInterpolation_dUdz(particlePosition, gridOrigin, gridSize, imax, jmax, kmax, dz, U);
    Real dVdx = trilinearInterpolation_dVdx(particlePosition, gridOrigin, gridSize, imax, jmax, kmax, dx, V);
    Real dVdz = trilinearInterpolation_dVdz(particlePosition, gridOrigin, gridSize, imax, jmax, kmax, dz, V);
    Real dWdx = trilinearInterpolation_dWdx(particlePosition, gridOrigin, gridSize, imax, jmax, kmax, dx, W);
    Real dWdy = trilinearInterpolation_dWdy(particlePosition, gridOrigin, gridSize, imax, jmax, kmax, dy, W);
    return rvec3(dWdy - dVdz, dUdz - dWdx, dVdx - dUdy);
}
