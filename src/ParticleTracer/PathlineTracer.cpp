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

#include "Intersection.hpp"
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

        // Break if the position is outside of the domain.
        if (glm::any(glm::lessThan(particlePosition, gridOrigin))
                || glm::any(glm::greaterThan(particlePosition, gridOrigin+gridSize))) {
            // Clamp the position to the boundary.
            rvec3 rayOrigin = particlePosition;
            rvec3 rayDirection = newParticlePosition - rayOrigin;
            float tNear, tFar;
            rayBoxIntersection(rayOrigin, rayDirection, gridOrigin, gridOrigin + gridSize, tNear, tFar);

            rvec3 boundaryParticlePosition = rayOrigin + tNear * rayDirection;
            trajectories.at(i).positions.push_back(boundaryParticlePosition);
            pushTrajectoryAttributes(
                    trajectories.at(i), gridOrigin, gridSize, imax, jmax, kmax, dx, dy, dz, U, V, W, P, T);
        } else {
            trajectories.at(i).positions.push_back(newParticlePosition);
            pushTrajectoryAttributes(
                    trajectories.at(i), gridOrigin, gridSize, imax, jmax, kmax, dx, dy, dz, U, V, W, P, T);
        }
    }
}

Trajectories PathlineTracer::getTrajectories(
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz,
        Real *U, Real *V, Real *W, Real *P, Real *T) {
    Trajectories trajectoriesCopy = trajectories;
    trajectories.clear();
    return trajectoriesCopy;
}
