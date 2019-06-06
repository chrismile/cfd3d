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

#ifndef CFD3D_FLAG_HPP
#define CFD3D_FLAG_HPP

#include <string>
#include <vector>
#include "Defines.hpp"

// Inline flag testing functions
// TODO: Add more inline functions.
inline bool isFluid(unsigned int flag) { return (flag >> 0) & 1; }


// Utility functions

/**
 * Initializes a flag array from a binary geometry file.
 * @param scenarioName The name of the CFD scenario to simulate.
 * @param geometryFilename The file name of the geometry file to load.
 * @param imax Number of cells in x direction inside of the domain.
 * @param jmax Number of cells in y direction inside of the domain.
 * @param kmax Number of cells in z direction inside of the domain.
 * @param FlagPtr A reference to a Flag array. This function uses new[] to allocate the memory of the Flag array.
 * The user is responsible for freeing the memory.
 */
void initFlagFromGeometryFile(const std::string &scenarioName, const std::string &geometryFilename,
        int imax, int jmax, int kmax, FlagType *&Flag);

/**
 * For scenarios not using arbitrary geometries.
 * This function justs sets the whole interior domain to fluid cells and all boundary cells to no-slip cells.
 * @param scenarioName The name of the CFD scenario to simulate.
 * @param imax Number of cells in x direction inside of the domain.
 * @param jmax Number of cells in y direction inside of the domain.
 * @param kmax Number of cells in z direction inside of the domain.
 * @param FlagPtr A reference to a Flag array. This function uses new[] to allocate the memory of the Flag array.
 */
void initFlagNoObstacles(const std::string &scenarioName, int imax, int jmax, int kmax, FlagType *&Flag);

/**
 * Upsamples (or downsamples if necessary) the passed values to match a certain resolution.
 * @param bitmapIn The bitmap to upscale.
 * @param widthIn The width of the bitmap to upscale.
 * @param heightIn The height of the bitmap to upscale.
 * @param bitmapOut The result.
 * @param widthOut The width of the resulting bitmap.
 * @param heightOut The height of the resulting bitmap.
 */
void nearestNeighborUpsampling(
        const std::vector<unsigned int> valuesIn,
        int widthIn,
        int heightIn,
        std::vector<unsigned int> &valuesOut,
        int widthOut,
        int heightOut);

#endif //CFD3D_FLAG_HPP
