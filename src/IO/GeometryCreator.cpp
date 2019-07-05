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

#include <boost/algorithm/string/predicate.hpp>
#include "PgmFile.hpp"
#include "GeometryFile.hpp"
#include "GeometryCreator.hpp"

/**
 * For more details regarding these values, please see docs/FlagsBitfield.pdf.
 */
enum GeometryValues {
    G_NO_SLIP, G_FREE_SLIP, G_OUTFLOW, G_INFLOW, G_FLUID,
    G_NO_SLIP_HOT, G_FREE_SLIP_HOT, G_INFLOW_HOT,
    G_NO_SLIP_COLD, G_FREE_SLIP_COLD, G_INFLOW_COLD,
    G_COUPLING
};

GeometryCreator::GeometryCreator(int imax, int jmax, int kmax, unsigned int boundaryType)
        : imax(imax), jmax(jmax), kmax(kmax) {
    // Initialize the domain with fluid cells in the interior, and no-slip cells on the boundary.
    geometryValues.resize((imax+2)*(jmax+2)*(kmax+2), boundaryType);

    #pragma omp parallel for
    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            for (int k = 1; k <= kmax; k++) {
                geometryValues[IDXFLAG(i,j,k)] = G_FLUID;
            }
        }
    }
}

void GeometryCreator::writeToBinGeoFile(const std::string &filename) {
    storeValuesInGeometryFile(filename, geometryValues, imax + 2, jmax + 2, kmax + 2);
}

void GeometryCreator::layersFromPgmFile(const std::string &filename, int layerStart, int layerEnd) {
    int pgmWidth, pgmHeight;
    std::vector<unsigned int> pgmValuesRead = loadPgmFile(filename, &pgmWidth, &pgmHeight);
    std::vector<unsigned int> pgmValues;
    nearestNeighborUpsampling2D(pgmValuesRead, pgmWidth, pgmHeight, pgmValues, imax+2, jmax+2);

    #pragma omp parallel for
    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            for (int k = layerStart; k <= layerEnd; k++) {
                geometryValues[IDXFLAG(i,j,k)] = pgmValues[i*(imax+2) + j];
            }
        }
    }
}

void GeometryCreator::setLayersConstant(FlagType value, int layerStart, int layerEnd) {
    #pragma omp parallel for
    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            for (int k = layerStart; k <= layerEnd; k++) {
                geometryValues[IDXFLAG(i,j,k)] = value;
            }
        }
    }
}

void GeometryCreator::setLayersInObject(
        FlagType value, int layerStart, int layerEnd, std::function<bool(int,int,int)> membershipFunctor) {
    #pragma omp parallel for
    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            for (int k = layerStart; k <= layerEnd; k++) {
                if (membershipFunctor(i,j,k)) {
                    geometryValues[IDXFLAG(i,j,k)] = value;
                }
            }
        }
    }
}


void createNaturalConvectionGeometry(
        const std::string &scenarioName, const std::string &geometryFilename, int imax, int jmax, int kmax) {
    GeometryCreator geometryCreator(imax, jmax, kmax, G_NO_SLIP);
    // Hot wall: Left.
    geometryCreator.setLayersInObject(G_NO_SLIP_HOT, 0, kmax+1, [&](int i, int j, int k) {
        return i == 0 && j > 0 && k > 0 && j <= jmax && k <= kmax;
    });
    // Cold wall: Right.
    geometryCreator.setLayersInObject(G_NO_SLIP_COLD, 0, kmax+1, [&](int i, int j, int k) {
        return i == imax+1 && j > 0 && k > 0 && j <= jmax && k <= kmax;
    });
    geometryCreator.writeToBinGeoFile(geometryFilename);
}

void createRayleighBenardGeometry(
        const std::string &scenarioName, const std::string &geometryFilename, int imax, int jmax, int kmax) {
    GeometryCreator geometryCreator(imax, jmax, kmax, G_NO_SLIP);
    // Hot wall: Bottom.
    geometryCreator.setLayersInObject(G_NO_SLIP_HOT, 0, kmax+1, [&](int i, int j, int k) {
        return j == 0 && i > 0 && k > 0 && i <= imax && k <= kmax;
    });
    // Cold wall: Top.
    geometryCreator.setLayersInObject(G_NO_SLIP_COLD, 0, kmax+1, [&](int i, int j, int k) {
        return j == jmax+1 && i > 0 && k > 0 && i <= imax && k <= kmax;
    });
    geometryCreator.writeToBinGeoFile(geometryFilename);
}

void createInflowTest(
        const std::string &scenarioName, const std::string &geometryFilename, int imax, int jmax, int kmax) {
    GeometryCreator geometryCreator(imax, jmax, kmax, G_FREE_SLIP);
    geometryCreator.setLayersInObject(G_INFLOW, 0, kmax+1, [&](int i, int j, int k) {
        return i == 0 && j >= 1 && k >= 1 && j <= jmax && k <= kmax;
    });
    geometryCreator.setLayersInObject(G_OUTFLOW, 0, kmax+1, [&](int i, int j, int k) {
        return i == imax+1 && j >= 1 && k >= 1 && j <= jmax && k <= kmax;
    });
    geometryCreator.writeToBinGeoFile(geometryFilename);
}

void generateScenario(
        const std::string &scenarioName, const std::string &geometryFilename, int imax, int jmax, int kmax) {
    if (scenarioName == "natural_convection") {
        createNaturalConvectionGeometry(scenarioName, geometryFilename, imax, jmax, kmax);
    } else if (boost::starts_with(scenarioName, "rayleigh_benard")) {
        createRayleighBenardGeometry(scenarioName, geometryFilename, imax, jmax, kmax);
    } else if (boost::starts_with(scenarioName, "inflow_test")) {
        createInflowTest(scenarioName, geometryFilename, imax, jmax, kmax);
    }
}