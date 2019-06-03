//
// Created by christoph on 03.06.19.
//

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
 * @param geometryFilename The filename of the .bingeo file to load.
 * @param width The width of the domain (including the boundary cells).
 * @param width The width of the domain (including the boundary cells).
 * @param height The height of the domain (including the boundary cells).
 * @param depth The depth of the domain (including the boundary cells).
 * @return The unsigned integer values stored in the file.
 */
void storeValuesInGeometryFile(
        const std::string &geometryFilename, const std::vector<uint32_t> &values, int width, int height, int depth);

#endif //CFD3D_GEOMETRYFILE_HPP
