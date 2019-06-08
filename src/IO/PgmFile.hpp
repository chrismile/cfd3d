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

#ifndef CFD3D_PGMFILE_HPP
#define CFD3D_PGMFILE_HPP

#include <string>
#include <vector>

/**
 * Loads a .pgm (portable gray map) image file.
 * @param filename The name of the file.
 * @param width A pointer to the field where the width of the image should be stored.
 * @param height A pointer to the field where the height of the image should be stored.
 * @return The grayscale values of the image stored as an array of 32-bit unsigned integer values.
 */
std::vector<unsigned int> loadPgmFile(const std::string &filename, int *width, int *height);

/**
 * Upsamples (or downsamples if necessary) a passed bitmap to match a certain resolution.
 * @param bitmapIn The bitmap to upscale.
 * @param widthIn The width of the bitmap to upscale.
 * @param heightIn The height of the bitmap to upscale.
 * @param bitmapOut The result.
 * @param widthOut The width of the resulting bitmap.
 * @param heightOut The height of the resulting bitmap.
 */
void nearestNeighborUpsampling2D(
        const std::vector<unsigned int> bitmapIn,
        int widthIn,
        int heightIn,
        std::vector<unsigned int> &bitmapOut,
        int widthOut,
        int heightOut
);

#endif //CFD3D_PGMFILE_HPP
