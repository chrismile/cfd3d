/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2019, Christoph Neuhauser, Stefan Haas, Paul Ng
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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
