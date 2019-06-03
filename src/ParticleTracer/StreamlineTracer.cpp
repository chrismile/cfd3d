//
// Created by christoph on 02.06.19.
//

#include "StreamlineTracer.hpp"
#include "Defines.hpp"

void traceStreamlineParticle(
        std::vector<rvec3> &currentLinePoints, const rvec3 &particleStartPosition,
        int imax, int jmax, int kmax, Real *U, Real *V, Real *W, Real dt) {
    currentLinePoints.push_back(particleStartPosition);
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
        rvec3 k1 = dt * getVectorVelocityAt(particlePosition, imax, jmax, kmax, U, V, W);
        rvec3 k2 = dt * getVectorVelocityAt(particlePosition + k1/Real(2.0), imax, jmax, kmax, U, V, W);
        rvec3 k3 = dt * getVectorVelocityAt(particlePosition + k2/Real(2.0), imax, jmax, kmax, U, V, W);
        rvec3 k4 = dt * getVectorVelocityAt(particlePosition + k3, imax, jmax, kmax, U, V, W);
        particlePosition = particlePosition + k1/Real(6.0) + k2/Real(3.0) + k3/Real(3.0) + k4/Real(6.0);

        // Add line segment between last and new position.
        currentLinePoints.push_back(particlePosition);
    } while(true);
}

std::vector<std::vector<rvec3>> StreamlineTracer::trace(
        const std::vector<rvec3> &particleSeedingLocations,
        int imax, int jmax, int kmax, Real *U, Real *V, Real *W, Real dt) {
    std::vector<std::vector<rvec3>> lines;
    lines.resize(particleSeedingLocations.size());
    for (size_t i = 0; i < particleSeedingLocations.size(); i++) {
        const rvec3 &particleStartPosition = particleSeedingLocations.at(i);
        std::vector<rvec3> &currentLinePoints = lines.at(i);
        traceStreamlineParticle(currentLinePoints, particleStartPosition, imax, jmax, kmax, U, V, W, dt);
    }
    return lines;
}
