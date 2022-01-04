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
#include <cstdlib>
#include "StringUtils.hpp"
#include "IOUtils.hpp"

/*
 * Fix for compatibility with older versions of GCC, as used e.g. on Ubuntu 18.04:
 * https://askubuntu.com/questions/1256440/how-to-get-libstdc-with-c17-filesystem-headers-on-ubuntu-18-bionic
 */
#if __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#error Missing the <filesystem> header.
#endif

void prepareOutputDirectory(
        const std::string &outputDirectory, const std::string &outputFilename, const std::string &outputFormatEnding,
        const std::string &lineDirectory, const std::string &geometryDirectory, bool shallWriteOutput) {
    // Create the directories.
    if (!fs::exists(outputDirectory)) {
        if (!fs::create_directory(outputDirectory)) {
            std::cerr << "Output directory could not be created." << std::endl;
            exit(1);
        }
    }
    if (!fs::exists(lineDirectory)) {
        if (!fs::create_directory(lineDirectory)) {
            std::cerr << "Line directory could not be created." << std::endl;
            exit(1);
        }
    }
    if (!fs::exists(geometryDirectory)) {
        if (!fs::create_directory(geometryDirectory)) {
            std::cerr << "Geometry directory could not be created." << std::endl;
            exit(1);
        }
    }

    // Delete all previous output files for the selected scenario.
    const std::string outputFilenameDot = outputFilename + ".";
    if (shallWriteOutput) {
        fs::path dir(outputDirectory);
        fs::directory_iterator end;
        for (fs::directory_iterator it(dir); it != end; ++it) {
            std::string filename = it->path().string();
#ifdef WIN32
            for (std::string::iterator it = filename.begin(); it != filename.end(); ++it) {
                if (*it == '\\') *it = '/';
            }
#endif
            if (endsWith(filename, outputFormatEnding) && startsWith(filename, outputFilenameDot)) {
                fs::remove(it->path());
            }
        }
    }
}
