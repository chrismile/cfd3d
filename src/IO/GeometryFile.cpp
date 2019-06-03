//
// Created by christoph on 03.06.19.
//

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
