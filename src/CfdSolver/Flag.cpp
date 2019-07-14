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

#include <cstring>
#include "IO/GeometryFile.hpp"
#include "Flag.hpp"

void nearestNeighborUpsampling(
        const std::vector<uint32_t> valuesIn,
        int widthIn,
        int heightIn,
        int depthIn,
        std::vector<uint32_t> &valuesOut,
        int widthOut,
        int heightOut,
        int depthOut
) {
    assert(widthIn >= 3 && heightIn >= 3 && depthIn >= 3 && widthOut >= 3 && heightOut >= 3 && depthOut >= 3);
    size_t size = static_cast<size_t>(widthOut * heightOut * depthOut);
    valuesOut.resize(size);

    // Size already matches?
    if (widthIn == widthOut && heightIn == heightOut && depthIn == depthOut) {
        memcpy(&valuesOut.front(), &valuesIn.front(), sizeof(uint32_t) * size);
        return;
    }

    // Copy the corners.
    for (int x = 0; x <= 1; x++) {
        for (int y = 0; y <= 1; y++) {
            for (int z = 0; z <= 1; z++) {
                valuesOut[static_cast<size_t>(x*(widthOut-1)*heightOut*depthOut + y*(heightOut-1)*depthOut + z*(depthOut-1))]
                        = valuesIn[static_cast<size_t>(x*(widthIn-1)*heightIn*depthIn + y*(heightIn-1)*depthIn + z*(depthIn-1))];
            }
        }
    }

    // Upscale borders in the two xy planes.
    for (int x = 1; x < widthOut - 1; x++) {
        for (int y = 1; y < heightOut - 1; y++) {
            int readX = (x - 1) * (widthIn - 2) / (widthOut - 2) + 1;
            int readY = (y - 1) * (heightIn - 2) / (heightOut - 2) + 1;
            valuesOut[static_cast<size_t>(x*heightOut*depthOut + y*depthOut + 0)]
                    = valuesIn[static_cast<size_t>(readX*heightIn*depthIn + readY*depthIn + 0)];
            valuesOut[static_cast<size_t>(x*heightOut*depthOut + y*depthOut + (depthOut-1))]
                    = valuesIn[static_cast<size_t>(readX*heightIn*depthIn + readY*depthIn + (depthIn-1))];
        }
    }

    // Upscale borders in the two xz planes.
    for (int x = 1; x < widthOut - 1; x++) {
        for (int z = 1; z < depthOut - 1; z++) {
            int readX = (x - 1) * (widthIn - 2) / (widthOut - 2) + 1;
            int readZ = (z - 1) * (depthIn - 2) / (depthOut - 2) + 1;
            valuesOut[static_cast<size_t>(x*heightOut*depthOut + 0*depthOut + z)]
                    = valuesIn[static_cast<size_t>(readX*heightIn*depthIn + 0*depthIn + readZ)];
            valuesOut[static_cast<size_t>(x*heightOut*depthOut + (heightOut-1)*depthOut + z)]
                    = valuesIn[static_cast<size_t>(readX*heightIn*depthIn + (heightIn-1)*depthIn + readZ)];
        }
    }

    // Upscale borders in the two yz planes.
    for (int y = 1; y < heightOut - 1; y++) {
        for (int z = 1; z < depthOut - 1; z++) {
            int readY = (y-1) * (heightIn-2) / (heightOut-2) + 1;
            int readZ = (z-1) * (depthIn-2) / (depthOut-2) + 1;
            valuesOut[static_cast<size_t>(0*heightOut*depthOut + y*depthOut + z)]
                    = valuesIn[static_cast<size_t>(0*heightIn*depthIn + readY*depthIn + readZ)];
            valuesOut[static_cast<size_t>((widthOut-1)*heightOut*depthOut + y*depthOut + z)]
                    = valuesIn[static_cast<size_t>((widthIn-1)*heightIn*depthIn + readY*depthIn + readZ)];
        }
    }

    // Upscale the interior cells.
    for (int x = 1; x < widthOut - 1; x++) {
        for (int y = 1; y < heightOut - 1; y++) {
            for (int z = 1; z < depthOut - 1; z++) {
                int readX = (x - 1) * (widthIn - 2) / (widthOut - 2) + 1;
                int readY = (y - 1) * (heightIn - 2) / (heightOut - 2) + 1;
                int readZ = (z - 1) * (depthIn - 2) / (depthOut - 2) + 1;
                valuesOut[static_cast<size_t>(x*heightOut*depthOut + y*depthOut + z)]
                        = valuesIn[static_cast<size_t>(readX*heightIn*depthIn + readY*depthIn + readZ)];
            }
        }
    }
}

void initFlagFromGeometryFile(const std::string &scenarioName, const std::string &geometryFilename,
        int imax, int jmax, int kmax, FlagType *&Flag) {
    int width, height, depth;
    std::vector<uint32_t> geometryValuesRead = loadValuesFromGeometryFile(geometryFilename, width, height, depth);
    std::vector<uint32_t> geometryValues;
    if (imax+2 != width || jmax+2 != height || kmax+2 != depth) {
        nearestNeighborUpsampling(geometryValuesRead, width, height, depth, geometryValues, imax+2, jmax+2, kmax+2);
    } else {
        geometryValues = geometryValuesRead;
    }

    // For an explanation of bit values, please see docs/BitfieldFlags.pdf.
    const unsigned int FLAG_LOOKUP_TABLE[] = {
            // no-slip, free-slip, outflow, inflow, fluid
            0x2, 0x4, 0x8, 0x10, 0x1,
            // no-slip (hot), free-slip (hot), inflow (hot)
            0x802, 0x804, 0x810,
            // no-slip (cold), free-slip (cold), inflow (cold)
            0x1002, 0x1004, 0x1010,
            // coupling
            0x2002
    };
    const unsigned int NUM_PGM_VALUES = sizeof(FLAG_LOOKUP_TABLE) / sizeof(*FLAG_LOOKUP_TABLE);

    #pragma omp parallel for
    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            for (int k = 0; k <= kmax+1; k++) {
                unsigned int pgmValue = geometryValues[IDXFLAG(i,j,k)];
                assert(pgmValue < NUM_PGM_VALUES);
                unsigned int flagValue = FLAG_LOOKUP_TABLE[pgmValue];
                Flag[IDXFLAG(i,j,k)] = flagValue;
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            for (int k = 0; k <= kmax+1; k++) {
                if (!isFluid(Flag[IDXFLAG(i,j,k)])) {
                    // Set B_L bit
                    if (i > 0 && (Flag[IDXFLAG(i - 1, j, k)] & 0x1) == 1) {
                        Flag[IDXFLAG(i, j, k)] |= 0x20;
                    }
                    // Set B_R bit
                    if (i <= imax && (Flag[IDXFLAG(i + 1, j, k)] & 0x1) == 1) {
                        Flag[IDXFLAG(i, j, k)] |= 0x40;
                    }
                    // Set B_D bit
                    if (j > 0 && (Flag[IDXFLAG(i, j - 1, k)] & 0x1) == 1) {
                        Flag[IDXFLAG(i, j, k)] |= 0x80;
                    }
                    // Set B_U bit
                    if (j <= jmax && (Flag[IDXFLAG(i, j + 1, k)] & 0x1) == 1) {
                        Flag[IDXFLAG(i, j, k)] |= 0x100;
                    }
                    // Set B_B bit
                    if (k > 0 && (Flag[IDXFLAG(i, j, k - 1)] & 0x1) == 1) {
                        Flag[IDXFLAG(i, j, k)] |= 0x200;
                    }
                    // Set B_F bit
                    if (k <= kmax && (Flag[IDXFLAG(i, j, k + 1)] & 0x1) == 1) {
                        Flag[IDXFLAG(i, j, k)] |= 0x400;
                    }
                }

                //#ifdef DEBUG
                int numNeighborsFluid = 0;
                if (i > 0) numNeighborsFluid += Flag[IDXFLAG(i-1,j,k)] & 0x1;
                if (i <= imax) numNeighborsFluid += Flag[IDXFLAG(i+1,j,k)] & 0x1;
                if (j > 0) numNeighborsFluid += Flag[IDXFLAG(i,j-1,k)] & 0x1;
                if (j <= jmax) numNeighborsFluid += Flag[IDXFLAG(i,j+1,k)] & 0x1;
                if (k > 0) numNeighborsFluid += Flag[IDXFLAG(i,j,k-1)] & 0x1;
                if (k <= kmax) numNeighborsFluid += Flag[IDXFLAG(i,j,k+1)] & 0x1;

                // Rule out boundary cells having more than three fluid neighbours, i.e.:
                // boundary cell => at most three fluid neighbors
                // Please remember rule for implication: a => b == not a or b
                assert((Flag[IDXFLAG(i,j,k)] & 0x1) == 1 || numNeighborsFluid <= 3);

                // Fluid cell may not have temperature source or boundary flag set.
                // fluid cell => none of the bits above set
                assert((Flag[IDXFLAG(i,j,k)] & 0x1) == 0 || (Flag[IDXFLAG(i,j,k)] & 0x381E) == 0);

                // Boundary cells with two opposite fluid cells are excluded (forbidden boundary cells).
                // boundary cell => ((fluid left => not fluid right) and (fluid down => not fluid up)
                // and (fluid back => not fluid front)).
                bool fluidLeft = i > 0 && (Flag[IDXFLAG(i-1,j,k)] & 0x1) == 1;
                bool fluidRight = i <= imax && (Flag[IDXFLAG(i+1,j,k)] & 0x1) == 1;
                bool fluidDown = j > 0 && (Flag[IDXFLAG(i,j-1,k)] & 0x1) == 1;
                bool fluidUp = j <= jmax && (Flag[IDXFLAG(i,j+1,k)] & 0x1) == 1;
                bool fluidBack = j <= jmax && (Flag[IDXFLAG(i,j,k-1)] & 0x1) == 1;
                bool fluidFront = j <= jmax && (Flag[IDXFLAG(i,j,k+1)] & 0x1) == 1;
                assert((Flag[IDXFLAG(i,j,k)] & 0x1) == 1
                        || ((!fluidLeft || !fluidRight) && (!fluidDown || !fluidUp) && (!fluidBack || !fluidFront)));

                // Only one of the four boundary type bits may be set.
                int boundaryNum = (Flag[IDXFLAG(i,j,k)] >> 1) & 0xF;
                assert(boundaryNum == 0x0 || boundaryNum == 0x1 || boundaryNum == 0x2
                        || boundaryNum == 0x4 || boundaryNum == 0x8);

                // Cells in the interior of the domain may only be fluid or no-slip.
                // interior cell => cell is fluid or no-slip
                assert(!(i > 0 && j > 0 && k > 0 && i <= imax && j <= jmax && k <= kmax)
                        || ((Flag[IDXFLAG(i,j,k)] & 0x1) == 1 || boundaryNum == 1));
                //#endif
            }
        }
    }
}

void initFlagNoObstacles(const std::string &scenarioName, int imax, int jmax, int kmax, FlagType *&Flag) {
    // Standard geometry: Fluid surrounded by no-slip boundary cells.
    #pragma omp parallel for
    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            for (int k = 0; k <= kmax+1; k++) {
                if (i == 0 || j == 0 || k == 0 || i == imax+1 || j == jmax+1 || k == kmax+1) {
                    Flag[IDXFLAG(i,j,k)] = 0x2; // no-slip
                } else {
                    Flag[IDXFLAG(i,j,k)] = 0x1; // fluid
                }
            }
        }
    }

    // Set neighbor values.
    #pragma omp parallel for
    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            for (int k = 0; k <= kmax+1; k++) {
                if (!isFluid(Flag[IDXFLAG(i,j,k)])) {
                    // Set B_L bit
                    if (i > 0 && (Flag[IDXFLAG(i - 1, j, k)] & 0x1) == 1) {
                        Flag[IDXFLAG(i, j, k)] |= 0x20;
                    }
                    // Set B_R bit
                    if (i <= imax && (Flag[IDXFLAG(i + 1, j, k)] & 0x1) == 1) {
                        Flag[IDXFLAG(i, j, k)] |= 0x40;
                    }
                    // Set B_D bit
                    if (j > 0 && (Flag[IDXFLAG(i, j - 1, k)] & 0x1) == 1) {
                        Flag[IDXFLAG(i, j, k)] |= 0x80;
                    }
                    // Set B_U bit
                    if (j <= jmax && (Flag[IDXFLAG(i, j + 1, k)] & 0x1) == 1) {
                        Flag[IDXFLAG(i, j, k)] |= 0x100;
                    }
                    // Set B_B bit
                    if (k > 0 && (Flag[IDXFLAG(i, j, k - 1)] & 0x1) == 1) {
                        Flag[IDXFLAG(i, j, k)] |= 0x200;
                    }
                    // Set B_F bit
                    if (k <= kmax && (Flag[IDXFLAG(i, j, k + 1)] & 0x1) == 1) {
                        Flag[IDXFLAG(i, j, k)] |= 0x400;
                    }
                }
            }
        }
    }
}
