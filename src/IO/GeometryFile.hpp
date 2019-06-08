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

#ifndef CFD3D_GEOMETRYFILE_HPP
#define CFD3D_GEOMETRYFILE_HPP

#include <string>
#include <vector>
#include <cstdint>

/**
 * Reads the values from a binary geometry file (.bingeo) and returns them.
 * For a mapping from these values to Flag entries, @see docs/BitfieldFlags.pdf.
 * @param geometryFilename The filename of the .bingeo file to load.
 * @param width The width of the domain (including the boundary cells).
 * @param height The height of the domain (including the boundary cells).
 * @param depth The depth of the domain (including the boundary cells).
 * @return The unsigned integer values stored in the file.
 */
std::vector<uint32_t> loadValuesFromGeometryFile(
        const std::string &geometryFilename, int &width, int &height, int &depth);

/**
 * Stores the passed unsigned integer values in a binary geometry file (.bingeo).
 * For a mapping from these values to Flag entries, @see docs/BitfieldFlags.pdf.
 * @param geometryFilename The filename of the .bingeo file to write.
 * @param width The width of the domain (including the boundary cells).
 * @param width The width of the domain (including the boundary cells).
 * @param height The height of the domain (including the boundary cells).
 * @param depth The depth of the domain (including the boundary cells).
 */
void storeValuesInGeometryFile(
        const std::string &geometryFilename, const std::vector<uint32_t> &values, int width, int height, int depth);

#endif //CFD3D_GEOMETRYFILE_HPP
