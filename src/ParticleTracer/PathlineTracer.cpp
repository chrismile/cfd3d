//
// Created by christoph on 02.06.19.
//

#include "PathlineTracer.hpp"

void PathlineTracer::setParticleSeedingLocations(
        const rvec3 &gridOrigin, const rvec3 &gridSize, const std::vector<rvec3> &particleSeedingLocations) {
    this->particleSeedingLocations = particleSeedingLocations;
    this->gridOrigin = gridOrigin;
    this->gridSize = gridSize;
    numParticles = particleSeedingLocations.size();
}

void PathlineTracer::timeStep(
        Real t, Real dt, int imax, int jmax, int kmax, Real dx, Real dy, Real dz,
        Real *U, Real *V, Real *W, Real *P, Real *T) {
    // Initialize the trajectories first if this is the first time step.
    if (trajectories.size() == 0) {
        trajectories.resize(numParticles);
        for (size_t i = 0; i < numParticles; i++) {
            trajectories.at(i).positions.push_back(particleSeedingLocations.at(i));
            pushTrajectoryAttributes(
                    trajectories.at(i), gridOrigin, gridSize, imax, jmax, kmax, dx, dy, dz, U, V, W, P, T);
        }
    }

    for (size_t i = 0; i < numParticles; i++) {
        rvec3 particlePosition = trajectories.at(i).positions.back();
        rvec3 newParticlePosition = integrateParticlePositionEuler(
                particlePosition, gridOrigin, gridSize,
                imax, jmax, kmax, U, V, W, dt);
        trajectories.at(i).positions.push_back(newParticlePosition);
        pushTrajectoryAttributes(
                trajectories.at(i), gridOrigin, gridSize, imax, jmax, kmax, dx, dy, dz, U, V, W, P, T);
    }
}

Trajectories PathlineTracer::getTrajectories() {
    Trajectories trajectoriesCopy = trajectories;
    trajectories.clear();
    return trajectoriesCopy;
}
