//
// Created by christoph on 01.06.19.
//

#ifndef CFD3D_PARTICLETRACER_HPP
#define CFD3D_PARTICLETRACER_HPP

#include <vector>
#include <glm/glm.hpp>
#include "Defines.hpp"

/**
 * SteadyFlowParticleTracer is the super-class of StreamlineTracer.
 * TimeVaryingParticleTracer is the super-class of PathlineTracer, StreaklineTracer.
 * They are used for tracing the paths of particles through the fluid in order to create the characteristic lines
 * of our flow (sometimes also called field lines).
 *
 * Characteristic lines are tangential to the flow.
 * => Line tangent = vector field direction
 *
 * dx(t) / dt = v(x(t), t)
 * x(0) = x_0
 *
 * Characterization of flow: Unsteady/time-varying vs steady flow.
 *
 * Types of characteristic lines of a flow:
 * - Path lines: Follow massless particles through time and space.
 * - Streak lines. Continuously release particles into flow at fixed position and connect particles.
 * - Stream lines: Trajectory of massless particles at one time step.
 * If the flow is steady, all types of characteristic lines are equal.
 *
 * Other possibilities for visualization in a future implementation: e.g. streak surfaces.
 * They are created by seeding particles along a curve and connecting them to form surfaces.
 */


/**
 * A super class for particle tracers for steady flows. This class is implemented by @see StreamlineTracer.
 * A steady flow is not dependent on time, thus can be executed on a single snapshot of U, V, and W forming the velocity
 * vector field.
 */
class SteadyFlowParticleTracer
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
            int imax, int jmax, int kmax, Real *U, Real *V, Real *W, Real dt)=0;
};

/**
 * A super class for particle tracers for time varying flows. This class is implemented by @see PathlineTracer and @see
 * StreaklineTracer.
 */
class TimeVaryingParticleTracer
{
public:
    /**
     * Sets the seeding positions of the particles to trace during the (time-dependent) simulation.
     * @param particleSeedingLocations The seeding locations of the particles to trace in staggered grid space.
     */
    virtual void setParticleSeedingLocations(const std::vector<rvec3> &particleSeedingLocations)=0;

    /**
     * Integrates the position of all particles with the passed time step size.
     * @param t The current time in the simulation.
     * @param dt The time step to use for integrating the particle position.
     * @param imax Number of cells in x direction inside of the domain.
     * @param jmax Number of cells in y direction inside of the domain.
     * @param kmax Number of cells in z direction inside of the domain.
     * @param U The velocities in x direction.
     * @param V The velocities in y direction.
     * @param W The velocities in z direction.
     */
    virtual void timeStep(Real t, Real dt, int imax, int jmax, int kmax, Real *U, Real *V, Real *W)=0;

    /**
     * Get the characteristic lines generated by the calls @see timeStep.
     * @return The characteristic lines (an array containing the arrays of line points).
     */
    virtual std::vector<std::vector<rvec3>> getLines()=0;
};

/**
 * Integrates a particle position in the given velocity field using the passed time step size with the explicit Euler
 * scheme.
 * @param particleStartPosition The position of the particle in staggered grid coordinates.
 * @param imax Number of cells in x direction inside of the domain.
 * @param jmax Number of cells in y direction inside of the domain.
 * @param kmax Number of cells in z direction inside of the domain.
 * @param U The velocities in x direction.
 * @param V The velocities in y direction.
 * @param W The velocities in z direction.
 * @param dt The time step to use for integrating the particle position.
 * @return The integrated particle position at time t+dt.
 */
rvec3 integrateParticlePositionEuler(const rvec3 &particlePosition, int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real dt);

/**
 * For internal use by the particle tracers. Returns the 3D velocity vector at a point in the staggered grid obtained
 * by trilinearly interpolating in U, V and W.
 * @param particlePosition
 * @param imax
 * @param jmax
 * @param kmax
 * @param U
 * @param V
 * @param W
 * @return
 */
rvec3 getVectorVelocityAt(const rvec3 &particlePosition, int imax, int jmax, int kmax, Real *U, Real *V, Real *W);

#endif //CFD3D_PARTICLETRACER_HPP
