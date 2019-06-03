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
     * @param particleSeedingLocations The seeding locations of the particles to trace in staggered grid coordinates.
     * @param imax Number of cells in x direction inside of the domain.
     * @param jmax Number of cells in y direction inside of the domain.
     * @param kmax Number of cells in z direction inside of the domain.
     * @param U The velocities in x direction.
     * @param V The velocities in y direction.
     * @param W The velocities in z direction.
     * @param dt The time step to use for integrating the particle position.
     * @return The characteristic lines (an array containing the arrays of line points).
     */
    virtual std::vector<std::vector<rvec3>> trace(const std::vector<rvec3> &particleSeedingLocations,
            int imax, int jmax, int kmax, Real *U, Real *V, Real *W, Real dt);
};

#endif //CFD3D_STREAMLINETRACER_HPP
