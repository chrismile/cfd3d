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
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include "IOUtils.hpp"

void prepareOutputDirectory(
        const std::string &outputDirectory, const std::string &outputFilename, const std::string &outputFormatEnding,
        const std::string &lineDirectory, const std::string &geometryDirectory, bool shallWriteOutput) {
    // Create the directories.
    if (!boost::filesystem::exists(outputDirectory)) {
        if (!boost::filesystem::create_directory(outputDirectory)) {
            std::cerr << "Output directory could not be created." << std::endl;
            exit(1);
        }
    }
    if (!boost::filesystem::exists(lineDirectory)) {
        if (!boost::filesystem::create_directory(lineDirectory)) {
            std::cerr << "Line directory could not be created." << std::endl;
            exit(1);
        }
    }
    if (!boost::filesystem::exists(geometryDirectory)) {
        if (!boost::filesystem::create_directory(geometryDirectory)) {
            std::cerr << "Geometry directory could not be created." << std::endl;
            exit(1);
        }
    }

    // Delete all previous output files for the selected scenario.
    if (shallWriteOutput) {
        boost::filesystem::path dir(outputDirectory);
        boost::filesystem::directory_iterator end;
        for (boost::filesystem::directory_iterator it(dir); it != end; ++it) {
            std::string filename = it->path().string();
#ifdef WIN32
            for (std::string::iterator it = filename.begin(); it != filename.end(); ++it) {
                if (*it == '\\') *it = '/';
            }
#endif
            if (boost::ends_with(filename, outputFormatEnding) && boost::starts_with(filename, outputFilename)) {
                boost::filesystem::remove(it->path());
            }
        }
    }
}