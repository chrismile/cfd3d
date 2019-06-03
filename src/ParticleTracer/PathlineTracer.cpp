//
// Created by christoph on 02.06.19.
//

#include "PathlineTracer.hpp"

void PathlineTracer::setParticleSeedingLocations(const std::vector<rvec3> &particleSeedingLocations) {
    this->particleSeedingLocations = particleSeedingLocations;
    numParticles = particleSeedingLocations.size();
    lines.resize(numParticles);
    for (size_t i = 0; i < numParticles; i++) {
        lines.at(i).push_back(particleSeedingLocations.at(i));
    }
}

void PathlineTracer::timeStep(Real t, Real dt, int imax, int jmax, int kmax, Real *U, Real *V, Real *W) {
    for (size_t i = 0; i < numParticles; i++) {
        rvec3 particlePosition = lines.at(i).back();
        rvec3 newParticlePosition = integrateParticlePositionEuler(particlePosition, imax, jmax, kmax, U, V, W, dt);
        lines.at(i).push_back(newParticlePosition);
    }
}

std::vector<std::vector<rvec3>> PathlineTracer::getLines() {
    std::vector<std::vector<rvec3>> linesCopy = lines;
    lines.clear();
    return linesCopy;
}
