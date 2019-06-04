//
// Created by christoph on 02.06.19.
//

#ifndef CFD3D_STREAMLINETRACER_HPP
#define CFD3D_STREAMLINETRACER_HPP

#include "ParticleTracer.hpp"

class StreamlineTracer : public SteadyFlowParticleTracer
{
public:
    /**
     * Traces the characteristic lines of a given steady velocity vector field.
     * @param particleSeedingLocations The seeding locations of the particles to trace in world space.
     * @param gridOrigin The origin of the grid in world coordinates.
     * @param gridSize The size of the grid (i.e. the extent in x, y and z) of the grid.
     * @param dt The time step to use for integrating the particle position.
     * @param imax Number of cells in x direction inside of the domain.
     * @param jmax Number of cells in y direction inside of the domain.
     * @param kmax Number of cells in z direction inside of the domain.
     * @param dx The cell size in x direction.
     * @param dy The cell size in y direction.
     * @param dz The cell size in z direction.
     * @param U The velocities in x direction.
     * @param V The velocities in y direction.
     * @param W The velocities in z direction.
     * @param P The pressure values.
     * @param T The temperature values.
     * @return The characteristic lines (an array containing the arrays of line points).
     */
    virtual Trajectories trace(
            const std::vector<rvec3> &particleSeedingLocations, const rvec3 &gridOrigin, const rvec3 &gridSize, Real dt,
            int imax, int jmax, int kmax, Real dx, Real dy, Real dz, Real *U, Real *V, Real *W, Real *P, Real *T);
};

#endif //CFD3D_STREAMLINETRACER_HPP
