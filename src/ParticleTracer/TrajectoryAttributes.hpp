//
// Created by christoph on 05.06.19.
//

#ifndef CFD3D_TRAJECTORYATTRIBUTES_HPP
#define CFD3D_TRAJECTORYATTRIBUTES_HPP

#include <vector>
#include <glm/glm.hpp>
#include "Defines.hpp"

struct Trajectory {
    std::vector<rvec3> positions;
    std::vector<std::vector<Real>> attributes;
};

typedef std::vector<Trajectory> Trajectories;

void pushTrajectoryAttributes(
        Trajectory &currentTrajectory, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz, Real *U, Real *V, Real *W, Real *P, Real *T);


/**
 * For internal use by the particle tracers. Returns the 3D velocity vector at a point in the staggered grid obtained
 * by trilinearly interpolating in U, V and W.
 * @param particleStartPosition The position of the particle in world coordinates.
 * @param gridOrigin The origin of the grid in world coordinates.
 * @param gridSize The size of the grid (i.e. the extent in x, y and z) of the grid.
 * @param imax Number of cells in x direction inside of the domain.
 * @param jmax Number of cells in y direction inside of the domain.
 * @param kmax Number of cells in z direction inside of the domain.
 * @param U The velocities in x direction.
 * @param V The velocities in y direction.
 * @param W The velocities in z direction.
 * @param dt The time step to use for integrating the particle position.
 * @return The velocity at the world space position particleStartPosition in the staggered grid.
 */
rvec3 getVectorVelocityAt(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real *U, Real *V, Real *W);

/**
 * curl v(p, t) = cross(nabla, v(p, t)) =
 * (d/dy v_z(p, t) - d/dz v_y(p, t),
 * d/dz v_x(p, t) - d/dx v_z(p, t),
 * d/dx v_y(p, t) - d/dy v_x(p, t))^T
 */
rvec3 getCurlAt(
        const rvec3 &particlePosition, const rvec3 &gridOrigin, const rvec3 &gridSize,
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz, Real *U, Real *V, Real *W);

#endif //CFD3D_TRAJECTORYATTRIBUTES_HPP
