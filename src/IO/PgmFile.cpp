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

#include <iostream>
#include <cassert>
#include <cstring>
#include "PgmFile.hpp"

std::vector<unsigned int> loadPgmFile(const std::string &filename, int &width, int &height, int &levels) {
    FILE *file = nullptr;
    if ((file = fopen(filename.c_str(), "rb")) == nullptr) {
        std::cerr << "Cannot open file \"" << filename << "\"." << std::endl;
        return {};
    }

    char line[1024];
    std::vector<unsigned int> pgmValues;

    if (fread(line, 1, 3, file) != 3) {
        fclose(file);
        std::cerr << "Wrong magic field in file \"" << filename << "\"." << std::endl;
        return {};
    }

    // Skip comment lines
    do {
        assert(fgets(line, sizeof(line), file) != nullptr);
    } while (*line == '#');

    // Read the image width and height
    sscanf(line, "%d %d\n", &width, &height);

    // Read the number of gray levels
    assert(fgets(line, sizeof line, file) != nullptr);
    sscanf(line, "%d\n", &levels);
    pgmValues.resize(size_t(width) * size_t(height));

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; ++x) {
            unsigned int value;
            if (fscanf(file, "%u", &value) == EOF) {
                fclose(file);
                std::cerr << "Invalid number of lines in file \"" << filename << "\"." << std::endl;
            }
            pgmValues[size_t(x) * size_t(height) + size_t(y)] = value;
        }
    }

    fclose(file);
    return pgmValues;
}

void nearestNeighborUpsamplingPgm2D(
        const std::vector<unsigned int> &bitmapIn,
        int widthIn,
        int heightIn,
        std::vector<unsigned int> &bitmapOut,
        int widthOut,
        int heightOut
) {
    assert(widthIn >= 3 && heightIn >= 3 && widthOut >= 3 && heightOut >= 3);
    size_t size = size_t(widthOut) * size_t(heightOut);
    bitmapOut.resize(size);

    // Size already matches?
    if (widthIn == widthOut && heightIn == heightOut) {
        memcpy(&bitmapOut.front(), &bitmapIn.front(), sizeof(unsigned int) * size);
        return;
    }

    // Copy corners
    bitmapOut[0*heightOut + 0] = bitmapIn[0*heightIn + 0];
    bitmapOut[(widthOut-1)*heightOut + 0] = bitmapIn[(widthIn-1)*heightIn + 0];
    bitmapOut[0*heightOut + heightOut-1] = bitmapIn[0 + heightIn-1];
    bitmapOut[(widthOut-1)*heightOut + heightOut-1] = bitmapIn[(widthIn-1)*heightIn + heightIn-1];

    // Upscale borders (in y direction)
    for (int y = 1; y < heightOut-1; y++) {
        int readY = (y-1) * (heightIn-2) / (heightOut-2) + 1;
        bitmapOut[0*heightOut + y] = bitmapIn[0*heightIn + readY];
        bitmapOut[(widthOut-1)*heightOut + y] = bitmapIn[(widthIn-1)*heightIn + readY];
    }

    // Upscale borders (in x direction)
    for (int x = 1; x < widthOut-1; x++) {
        int readX = (x-1) * (widthIn-2) / (widthOut-2) + 1;
        bitmapOut[x*heightOut + 0] = bitmapIn[readX*heightIn + 0];
        bitmapOut[x*heightOut + heightOut-1] = bitmapIn[readX*heightIn + heightIn-1];
    }

    // Upscale interior cells
    for (int y = 1; y < heightOut-1; y++) {
        for (int x = 1; x < widthOut-1; x++) {
            int readX = (x-1) * (widthIn-2) / (widthOut-2) + 1;
            int readY = (y-1) * (heightIn-2) / (heightOut-2) + 1;
            bitmapOut[x*heightOut + y*widthOut] = bitmapIn[readX*heightIn + readY];
        }
    }
}

void nearestNeighborUpsampling2D(
        const std::vector<unsigned int> &bitmapIn,
        int widthIn,
        int heightIn,
        std::vector<unsigned int> &bitmapOut,
        int widthOut,
        int heightOut
) {
    assert(widthIn >= 1 && heightIn >= 1 && widthOut >= 1 && heightOut >= 1);
    size_t size = size_t(widthOut) * size_t(heightOut);
    bitmapOut.resize(size);

    // Size already matches?
    if (widthIn == widthOut && heightIn == heightOut) {
        memcpy(&bitmapOut.front(), &bitmapIn.front(), sizeof(unsigned int) * size);
        return;
    }

    // Upscale interior cells
    for (int y = 0; y < heightOut; y++) {
        for (int x = 0; x < widthOut; x++) {
            int readX = x * widthIn / widthOut;
            int readY = y * heightIn / heightOut;
            bitmapOut[x*heightOut + y*widthOut] = bitmapIn[readX*heightIn + readY];
        }
    }
}
