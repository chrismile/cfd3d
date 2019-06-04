//
// Created by christoph on 02.06.19.
//

#include <algorithm>
#include "StreaklineTracer.hpp"

void StreaklineTracer::setParticleSeedingLocations(
        const rvec3 &gridOrigin, const rvec3 &gridSize, const std::vector<rvec3> &particleSeedingLocations) {
    this->particleSeedingLocations = particleSeedingLocations;
    this->gridOrigin = gridOrigin;
    this->gridSize = gridSize;
    lastInjectTime = 0.0;
    stepIndex = 0;
    numParticles = particleSeedingLocations.size();
}

void StreaklineTracer::timeStep(
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
        for (size_t j = 0; j < stepIndex; j++) {
            rvec3 particlePosition = trajectories.at(i).positions.at(j);
            trajectories.at(i).positions.at(j) = integrateParticlePositionEuler(
                    particlePosition, gridOrigin, gridSize,
                    imax, jmax, kmax, U, V, W, dt);
            pushTrajectoryAttributes(
                    trajectories.at(i), gridOrigin, gridSize, imax, jmax, kmax, dx, dy, dz, U, V, W, P, T);
        }
        if (t - lastInjectTime >= dtInject) {
            trajectories.at(i).positions.push_back(particleSeedingLocations.at(i));
            lastInjectTime = t;
        }
    }
    stepIndex++;
}

Trajectories StreaklineTracer::getTrajectories() {
    for (size_t i = 0; i < numParticles; i++) {
        std::reverse(trajectories.at(i).positions.begin(), trajectories.at(i).positions.end());
        std::reverse(trajectories.at(i).attributes.begin(), trajectories.at(i).attributes.end());
    }
    Trajectories trajectoriesCopy = trajectories;
    trajectories.clear();
    return trajectoriesCopy;
}
