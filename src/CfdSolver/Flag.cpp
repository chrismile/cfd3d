//
// Created by christoph on 01.06.19.
//

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
    for (int x = 1; x < heightOut - 1; x++) {
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
    for (int x = 1; x < heightOut - 1; x++) {
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
        int imax, int jmax, int kmax, FlagType *&FlagPtr) {
    int width, height, depth;
    std::vector<uint32_t> geometryValuesRead = loadValuesFromGeometryFile(geometryFilename, width, height, depth);
    std::vector<uint32_t> geometryValues;
    if (imax != width || jmax != height || kmax != depth) {
        nearestNeighborUpsampling(geometryValuesRead, width, height, depth, geometryValues, imax+2, jmax+2, kmax+2);
    }

    // For an explanation of bit values, please see docs/BitfieldFlags.pdf.
    const unsigned int FLAG_LOOKUP_TABLE[] = {
            // no-slip, free-slip, outflow, inflow, fluid
            0x2, 0x4, 0x8, 0x10, 0x1,
            // no-slip (hot), free-slip (hot), inflow (hot)
            0x202, 0x204, 0x210,
            // no-slip (cold), free-slip (cold), inflow (cold)
            0x402, 0x404, 0x410,
            // coupling
            0x802
    };

    // TODO: Implement.
}

void initFlagNoObstacles(const std::string &scenarioName, int imax, int jmax, int kmax, FlagType *&FlagPtr) {
    // TODO: Implement.
}
