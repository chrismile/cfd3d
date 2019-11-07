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
#include "Intersection.hpp"

void traceStreamlineParticle(
        Trajectory &currentTrajectory, const rvec3 &particleStartPosition, const rvec3 &gridOrigin,
        const rvec3 &gridSize, Real dt, int imax, int jmax, int kmax, Real dx, Real dy, Real dz,
        Real *U, Real *V, Real *W, Real *P, Real *T) {
    rvec3 particlePosition = particleStartPosition;
    rvec3 oldParticlePosition;
    glm::ivec3 particleGridPosition;

    int iterationCounter = 0;
    const int MAX_ITERATIONS = 2000;
    Real lineLength = 0.0;
    const Real MAX_LINE_LENGTH = gridSize.x + gridSize.y + gridSize.z;
    while (iterationCounter <= MAX_ITERATIONS && lineLength <= MAX_LINE_LENGTH) {
        oldParticlePosition = particlePosition;
        // Break if the position is outside of the domain.
        if (glm::any(glm::lessThan(particlePosition, gridOrigin))
                || glm::any(glm::greaterThan(particlePosition, gridOrigin+gridSize))) {
            if (currentTrajectory.positions.size() >= 1) {
                // Clamp the position to the boundary.
                rvec3 rayOrigin = currentTrajectory.positions.back();
                rvec3 rayDirection = particlePosition - rayOrigin;
                Real tNear, tFar;
                rayBoxIntersection(rayOrigin, rayDirection, gridOrigin, gridOrigin + gridSize, tNear, tFar);

                rvec3 boundaryParticlePosition = rayOrigin + tNear * rayDirection;
                currentTrajectory.positions.push_back(boundaryParticlePosition);
                pushTrajectoryAttributes(
                        currentTrajectory, gridOrigin, gridSize, imax, jmax, kmax, dx, dy, dz, U, V, W, P, T);
            }
            break;
        }

        // Add line segment between last and new position.
        currentTrajectory.positions.push_back(particlePosition);
        pushTrajectoryAttributes(
                currentTrajectory, gridOrigin, gridSize, imax, jmax, kmax, dx, dy, dz, U, V, W, P, T);

        // Integrate to the new position using Runge-Kutta of 4th order.
        /*
        rvec3 k1 = dt * getVectorVelocityAt(particlePosition, gridOrigin, gridSize, imax, jmax, kmax, U, V, W);
        rvec3 k2 = dt * getVectorVelocityAt(particlePosition + k1/Real(2.0), gridOrigin, gridSize, imax, jmax, kmax,
                U, V, W);
        rvec3 k3 = dt * getVectorVelocityAt(particlePosition + k2/Real(2.0), gridOrigin, gridSize, imax, jmax, kmax,
                U, V, W);
        rvec3 k4 = dt * getVectorVelocityAt(particlePosition + k3, gridOrigin, gridSize, imax, jmax, kmax, U, V, W);
        particlePosition = particlePosition + k1/Real(6.0) + k2/Real(3.0) + k3/Real(3.0) + k4/Real(6.0);
        */

        // Integrate to the new position using the Runge-Kutta-Fehlberg method (RKF45).
        // For more details see: http://maths.cnam.fr/IMG/pdf/RungeKuttaFehlbergProof.pdf
        rvec3 approximationRK4, approximationRK5;
        bool timestepNeedsAdaptation = false;
        do {
            rvec3 k1 = dt * getVectorVelocityAt(particlePosition, gridOrigin, gridSize, imax, jmax, kmax, U, V, W);
            rvec3 k2 = dt * getVectorVelocityAt(
                    particlePosition + k1*Real(1.0/4.0),
                    gridOrigin, gridSize, imax, jmax, kmax, U, V, W);
            rvec3 k3 = dt * getVectorVelocityAt(
                    particlePosition + k1*Real(3.0/32.0) + k2*Real(9.0/32.0),
                    gridOrigin, gridSize, imax, jmax, kmax, U, V, W);
            rvec3 k4 = dt * getVectorVelocityAt(
                    particlePosition + k1*Real(1932.0/2197.0) - k2*Real(7200.0/2197.0) + k3*Real(7296.0/2197.0),
                    gridOrigin, gridSize, imax, jmax, kmax, U, V, W);
            rvec3 k5 = dt * getVectorVelocityAt(
                    particlePosition + k1*Real(439.0/216.0) - k2*Real(8.0) + k3*Real(3680.0/513.0)
                    - k4*Real(845.0/4104.0),
                    gridOrigin, gridSize, imax, jmax, kmax, U, V, W);
            rvec3 k6 = dt * getVectorVelocityAt(
                    particlePosition - k1*Real(8.0/27.0) + k2*Real(2.0) - k3*Real(3544.0/2565.0)
                    + k4*Real(1859.0/4104.0) - k5*Real(11.0/40.0),
                    gridOrigin, gridSize, imax, jmax, kmax, U, V, W);
            approximationRK4 = particlePosition + k1*Real(25.0/216.0) + k3*Real(1408.0/2565.0)
                                     + k4*Real(2197.0/4101.0) - k5*Real(1.0/5.0);
            approximationRK5 = particlePosition + k1*Real(16.0/135.0) + k3*Real(6656.0/12825.0)
                                     + k4*Real(28561.0/56430.0) - k5*Real(9.0/50.0) + k6*Real(2.0/55.0);
            Real s = std::pow(tol*dt / (Real(2.0) * glm::length(approximationRK5 - approximationRK4)), Real(1.0/4.0));
            if (s < 0.9 || s > 1.1) {
                timestepNeedsAdaptation = true;
            }
            dt = s*dt;
        } while(timestepNeedsAdaptation);
        particlePosition = approximationRK4;


        Real segmentLength = glm::length(particlePosition - oldParticlePosition);
        lineLength += segmentLength;

        // Have we reached a singular point?
        if (segmentLength < Real(0.000001)) {
            break;
        }

        oldParticlePosition = particlePosition;
        iterationCounter++;
    }
}

Trajectories StreamlineTracer::trace(
        const std::vector<rvec3> &particleSeedingLocations, const rvec3 &gridOrigin, const rvec3 &gridSize, Real tolDt,
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz, Real *U, Real *V, Real *W, Real *P, Real *T) {
    Trajectories trajectories;
    trajectories.resize(particleSeedingLocations.size());
    for (size_t i = 0; i < particleSeedingLocations.size(); i++) {
        const rvec3 &particleStartPosition = particleSeedingLocations.at(i);
        Trajectory &currentTrajectory = trajectories.at(i);
        traceStreamlineParticle(
                currentTrajectory, particleStartPosition, gridOrigin, gridSize,
                dt, imax, jmax, kmax, dx, dy, dz, U, V, W, P, T);
    }
    return trajectories;
}
