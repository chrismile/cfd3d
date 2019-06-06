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

#include "StreamlineTracer.hpp"
#include "Defines.hpp"

void traceStreamlineParticle(
        Trajectory &currentTrajectory, const rvec3 &particleStartPosition, const rvec3 &gridOrigin,
        const rvec3 &gridSize, Real dt, int imax, int jmax, int kmax, Real dx, Real dy, Real dz,
        Real *U, Real *V, Real *W, Real *P, Real *T) {
    currentTrajectory.positions.push_back(particleStartPosition);
    rvec3 particlePosition = particleStartPosition;
    glm::ivec3 particleGridPosition;

    do {
        particleGridPosition = glm::ivec3(particlePosition);
        // Break if the position is outside of the domain.
        if (glm::any(glm::lessThan(particleGridPosition, glm::ivec3(1)))
            || glm::any(glm::greaterThan(particleGridPosition, glm::ivec3(imax, jmax, kmax)))) {
            break;
        }

        // Integrate to the new position using Runge-Kutta of 4th order.
        rvec3 k1 = dt * getVectorVelocityAt(particlePosition, gridOrigin, gridSize, imax, jmax, kmax, U, V, W);
        rvec3 k2 = dt * getVectorVelocityAt(particlePosition + k1/Real(2.0), gridOrigin, gridSize, imax, jmax, kmax,
                U, V, W);
        rvec3 k3 = dt * getVectorVelocityAt(particlePosition + k2/Real(2.0), gridOrigin, gridSize, imax, jmax, kmax,
                U, V, W);
        rvec3 k4 = dt * getVectorVelocityAt(particlePosition + k3, gridOrigin, gridSize, imax, jmax, kmax, U, V, W);
        particlePosition = particlePosition + k1/Real(6.0) + k2/Real(3.0) + k3/Real(3.0) + k4/Real(6.0);

        // Add line segment between last and new position.
        currentTrajectory.positions.push_back(particlePosition);
        pushTrajectoryAttributes(
                currentTrajectory, gridOrigin, gridSize, imax, jmax, kmax, dx, dy, dz, U, V, W, P, T);
    } while(true);
}

Trajectories StreamlineTracer::trace(
        const std::vector<rvec3> &particleSeedingLocations, const rvec3 &gridOrigin, const rvec3 &gridSize, Real dt,
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz, Real *U, Real *V, Real *W, Real *P, Real *T) {
    Trajectories trajectories;
    trajectories.resize(particleSeedingLocations.size());
    for (size_t i = 0; i < particleSeedingLocations.size(); i++) {
        const rvec3 &particleStartPosition = particleSeedingLocations.at(i);
        Trajectory &currentTrajectory = trajectories.at(i);
        traceStreamlineParticle(
                currentTrajectory, gridOrigin, gridOrigin, particleStartPosition,
                dt, imax, jmax, kmax, dx, dy, dz, U, V, W, P, T);
    }
    return trajectories;
}
