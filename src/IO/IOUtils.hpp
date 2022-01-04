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

#ifndef CFD3D_IOUTILS_HPP
#define CFD3D_IOUTILS_HPP

#include <string>

/**
 * Makes sure the directories 'outputDirectory', 'outputDirectory' and 'outputDirectory' exists.
 * @param outputDirectory The output directory name.
 * @param outputFilename The filename of the output file (not including the file ending).
 * @param outputFormatEnding The file ending of the format to use.
 * @param lineDirectory The line directory name.
 * @param geometryDirectory The geometry directory name.
 */
void prepareOutputDirectory(
        const std::string &outputDirectory, const std::string &outputFilename, const std::string &outputFormatEnding,
        const std::string &lineDirectory, const std::string &geometryDirectory, bool shallWriteOutput);

/**
 * Returns whether the passed file exists.
 * @param filename The path to the file.
 * @return True if the file exists and false otherwise.
 */
bool fileExists(const std::string& filename);

#endif //CFD3D_IOUTILS_HPP
