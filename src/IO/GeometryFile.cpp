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

#include <fstream>
#include <iostream>
#include "BinaryStream.hpp"
#include "GeometryFile.hpp"

const uint32_t GEOMETRY_FILE_FORMAT_VERSION = 1u;

std::vector<uint32_t> loadValuesFromGeometryFile(
        const std::string &geometryFilename, int &width, int &height, int &depth) {
    std::ifstream file(geometryFilename.c_str(), std::ifstream::binary);
    if (!file.is_open()) {
        std::cerr << "Error in loadValuesFromGeometryFile: File \"" + geometryFilename + "\" not found." << std::endl;
        return {};
    }

    file.seekg(0, file.end);
    size_t size = file.tellg();
    file.seekg(0);
    char *buffer = new char[size]; // BinaryReadStream does deallocation
    file.read(buffer, size);
    file.close();

    BinaryReadStream stream(buffer, size);
    uint32_t version;
    stream.read(version);
    if (version != GEOMETRY_FILE_FORMAT_VERSION) {
        std::cerr << "Error in readMesh3D: Invalid version in file \"" << geometryFilename << "\"." << std::endl;
        return {};
    }

    uint32_t uwidth, uheight, udepth;
    stream.read(uwidth);
    stream.read(uheight);
    stream.read(udepth);
    width = static_cast<int>(uwidth);
    height = static_cast<int>(uheight);
    depth = static_cast<int>(udepth);

    std::vector<uint32_t> values;
    values.resize(uwidth*uheight*udepth);
    stream.read((void*)&values.front(), sizeof(uint32_t)*values.size());
    return values;
}

void storeValuesInGeometryFile(
        const std::string &geometryFilename, const std::vector<uint32_t> &values, int width, int height, int depth) {
    std::ofstream file(geometryFilename.c_str(), std::ofstream::binary);
    if (!file.is_open()) {
        std::cerr << "Error in storeValuesInGeometryFile: File \"" + geometryFilename + "\" not found." << std::endl;
        return;
    }
    assert(values.size() == width*height*depth);

    BinaryWriteStream stream;
    stream.write((uint32_t)GEOMETRY_FILE_FORMAT_VERSION);
    stream.write((uint32_t)width);
    stream.write((uint32_t)height);
    stream.write((uint32_t)depth);
    stream.write((void*)&values.front(), sizeof(uint32_t)*values.size());

    file.write((const char*)stream.getBuffer(), stream.getSize());
    file.close();
}
