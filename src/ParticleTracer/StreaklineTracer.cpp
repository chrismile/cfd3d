//
// Created by christoph on 02.06.19.
//

#include <algorithm>
#include "StreaklineTracer.hpp"

void StreaklineTracer::setParticleSeedingLocations(const std::vector<rvec3> &particleSeedingLocations) {
    this->particleSeedingLocations = particleSeedingLocations;
    stepIndex = 0;
    numParticles = particleSeedingLocations.size();
    lines.resize(numParticles);
    for (size_t i = 0; i < numParticles; i++) {
        lines.at(i).push_back(particleSeedingLocations.at(i));
    }
}

void StreaklineTracer::timeStep(Real t, Real dt, int imax, int jmax, int kmax, Real *U, Real *V, Real *W) {
    for (size_t i = 0; i < numParticles; i++) {
        for (size_t j = 0; j < stepIndex; j++) {
            rvec3 particlePosition = lines.at(i).at(j);
            lines.at(i).at(j) = integrateParticlePositionEuler(particlePosition, imax, jmax, kmax, U, V, W, dt);
        }
        lines.at(i).push_back(particleSeedingLocations.at(i));
    }
    stepIndex++;
}

std::vector<std::vector<rvec3>> StreaklineTracer::getLines() {
    for (size_t i = 0; i < numParticles; i++) {
        std::reverse(lines.at(i).begin(), lines.at(i).end());
    }
    std::vector<std::vector<rvec3>> linesCopy = lines;
    lines.clear();
    return linesCopy;
}
