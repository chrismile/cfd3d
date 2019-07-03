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

#include <glm/glm.hpp>
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

GeometryCreator::GeometryCreator(int imax, int kmax, int jmax, unsigned int boundaryType)
        : imax(imax), jmax(jmax), kmax(kmax) {
    // Initialize the domain with fluid cells in the interior, and no-slip cells on the boundary.
    geometryValues.resize((imax+2)*(jmax+2)*(kmax+2), boundaryType);
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
    nearestNeighborUpsamplingPgm2D(pgmValuesRead, pgmWidth, pgmHeight, pgmValues, imax+2, jmax+2);

    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            for (int k = layerStart; k <= layerEnd; k++) {
                geometryValues[IDXFLAG(i,j,k)] = pgmValues[i*(imax+2) + j];
            }
        }
    }
}

void GeometryCreator::setLayersConstant(FlagType value, int layerStart, int layerEnd) {
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

bool isValidCell(int i, int j, int k, int imax, int jmax, int kmax, FlagType *Flag) {
    // Boundary cells with two opposite fluid cells are excluded (forbidden boundary cells).
    // boundary cell => ((fluid left => not fluid right) and (fluid bottom => not fluid top)
    // and (fluid back => not fluid front))
    bool fluidLeft = i > 0 && Flag[IDXFLAG(i-1,j,k)] == G_FLUID;
    bool fluidRight = i <= imax && Flag[IDXFLAG(i+1,j,k)] == G_FLUID;
    bool fluidBottom = j > 0 && Flag[IDXFLAG(i,j-1,k)] == G_FLUID;
    bool fluidTop = j <= jmax && Flag[IDXFLAG(i,j+1,k)] == G_FLUID;
    bool fluidBack = j > 0 && Flag[IDXFLAG(i,j,k-1)] == G_FLUID;
    bool fluidFront = j <= jmax && Flag[IDXFLAG(i,j,k+1)] == G_FLUID;
    return (Flag[IDXFLAG(i,j,k)] == G_FLUID
            || ((!fluidLeft || !fluidRight) && (!fluidBottom || !fluidTop)
                && (!fluidBack || !fluidFront)));
}

void GeometryCreator::removeInvalidCells() {
    // Go over all layers in y direction.
    for (int j = 1; j <= jmax+1; j++) {
        for (int i = 1; i <= imax+1; i++) {
            for (int k = 1; k <= kmax+1; k++) {
                if (!isValidCell(i, j, k, imax, jmax, kmax, &geometryValues.front())) {
                    geometryValues[IDXFLAG(i,j,k)] = G_FLUID;
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

void createFlowOverStepGeometry(
        const std::string &scenarioName, const std::string &geometryFilename, int imax, int jmax, int kmax) {
    GeometryCreator geometryCreator(imax, jmax, kmax, G_NO_SLIP);
    // Step box.
    geometryCreator.setLayersInObject(G_NO_SLIP, 0, kmax+1, [&](int i, int j, int k) {
        return i <= jmax/2 && j <= jmax/2;
    });
    geometryCreator.writeToBinGeoFile(geometryFilename);
}

void createTowerGeometry(
        const std::string &scenarioName, const std::string &geometryFilename, int imax, int jmax, int kmax) {
    // Create a tower centered at 'towerCenter' with the specified width (= depth) and height.
    glm::ivec3 towerCenter(imax/3, 0, kmax/2);
    int width = kmax / 4;
    int height = jmax*2/3;
    glm::ivec3 towerMin = towerCenter + glm::ivec3(-width/2, 0, -width/2);
    glm::ivec3 towerMax = towerCenter + glm::ivec3(width/2, height, width/2);

    GeometryCreator geometryCreator(imax, jmax, kmax, G_NO_SLIP);
    geometryCreator.setLayersInObject(G_NO_SLIP, 0, kmax+1, [&](int i, int j, int k) {
        glm::ivec3 coords(i, j, k);
        return glm::all(glm::greaterThanEqual(coords, towerMin)) && glm::all(glm::lessThanEqual(coords, towerMax));
    });
    geometryCreator.writeToBinGeoFile(geometryFilename);
}

void createTerrainGeometry(
        const std::string &scenarioName, const std::string &geometryFilename, int imax, int jmax, int kmax) {
    // Load the height map from the specified .pgm image file.
    int pgmWidth = 0;
    int pgmHeight = 0;
    std::vector<unsigned int> heightMapPgmData = loadPgmFile("../geometry-pgm/heightmap1.pgm", &pgmWidth, &pgmHeight);
    std::vector<unsigned int> heightMapInt;
    nearestNeighborUpsampling2D(heightMapPgmData, pgmWidth, pgmHeight, heightMapPgmData, imax, kmax);
    std::vector<float> heightMapFloat;
    heightMapFloat.resize(heightMapPgmData.size());
    for (size_t i = 0; i < heightMapPgmData.size(); i++) {
        heightMapFloat[i] = heightMapPgmData[i]/255.0f;
    }

    // Now set all cells below the specified heights to no-slip obstacles.
    GeometryCreator geometryCreator(imax, jmax, kmax, G_NO_SLIP);
    geometryCreator.setLayersInObject(G_NO_SLIP, 1, kmax, [&](int i, int j, int k) {
        if (i < 1 || i > imax || j < 1 || j > jmax) {
            return false;
        }
        float heightMapEntry = heightMapFloat.at(i * kmax + k);
        int integerHeight = static_cast<int>(heightMapEntry * jmax);
        return j <= integerHeight;
    });
    geometryCreator.removeInvalidCells();
    geometryCreator.writeToBinGeoFile(geometryFilename);
}

void generateScenario(
        const std::string &scenarioName, const std::string &geometryFilename, int imax, int jmax, int kmax) {
    if (scenarioName == "natural_convection") {
        createNaturalConvectionGeometry(scenarioName, geometryFilename, imax, jmax, kmax);
    } else if (boost::starts_with(scenarioName, "rayleigh_benard")) {
        createRayleighBenardGeometry(scenarioName, geometryFilename, imax, jmax, kmax);
    } else if (boost::starts_with(scenarioName, "flow_over_step")) {
        createFlowOverStepGeometry(scenarioName, geometryFilename, imax, jmax, kmax);
    } else if (boost::starts_with(scenarioName, "single_tower")) {
        createFlowOverStepGeometry(scenarioName, geometryFilename, imax, jmax, kmax);
    }
}