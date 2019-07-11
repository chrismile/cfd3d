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

#ifndef CFD3D_GEOMETRYCREATOR_HPP
#define CFD3D_GEOMETRYCREATOR_HPP

#include <string>
#include <functional>
#include "Defines.hpp"

/**
 * GeometryCreator is a class for programmatically generating .bingeo files storing the geometry information for
 * different scenarios. It offers functions to set certain layers in the domain, to set layers from two-dimensional
 * .pgm files, and to use object stencils for overwriting values in the domain.
 * When talking about 'layers', two-dimensional k-slices of our i,j,k domain are meant.
 */
class GeometryCreator
{
public:
    /**
     * Initializes the domain of the geometry data.
     * For the inside of the domain, fluid cells are used for initialization.
     * For the boundary, the passed boundary type is used.
     * @param imax Number of cells in x direction inside of the domain.
     * @param jmax Number of cells in y direction inside of the domain.
     * @param kmax Number of cells in z direction inside of the domain.
     * @param boundaryType The type of boundary to use for initialization (usually no-slip or free-slip).
     */
    GeometryCreator(int imax, int jmax, int kmax, unsigned int boundaryType);

    /**
     * Writes the geometry values to a binary geometry file (.bingeo). For more details @see storeValuesInGeometryFile.
     * @param filename The filename of the .bingeo file.
     */
    void writeToBinGeoFile(const std::string &filename);

    /**
     * Sets the values of one or multiple layers to the values stored in a .pgm file. For more details @see loadPgmFile.
     * If the size of the .pgm file does not match, nearest neighbor upsampling/downsampling is used for the interior
     * of the domain and the boundaries are copied to the larger resolution.
     * @param filename The filename of the .pgm file.
     * @param layerStart The first layer (inclusive) to overwrite with the .pgm values.
     * @param layerEnd The last layer (inclusive) to overwrite with the .pgm values.
     */
    void layersFromPgmFile(const std::string &filename, int layerStart, int layerEnd);

    /**
     * Sets certain layers to a constant value.
     * @param value The value to use for the specified layers.
     * @param layerStart The first layer (inclusive) to overwrite with the value.
     * @param layerEnd The last layer (inclusive) to overwrite with the value.
     */
    void setLayersConstant(FlagType value, int layerStart, int layerEnd);

    /**
     * Overwrites the geometry values at all points passing the membership functor.
     * @param value The value to use for overwriting.
     * @param layerStart The first layer (inclusive) to test.
     * @param layerEnd The last layer (inclusive) to test.
     * @param membershipFunctor If this function returns true for a position (i,j,k), the value at this position
     * gets overwritten with the passed value.
     */
    void setLayersInObject(
            FlagType value, int layerStart, int layerEnd, std::function<bool(int,int,int)> membershipFunctor);

    /**
     * Removes invalid obstacle cells and replaces them with fluid cells.
     */
    void removeInvalidCells();

    /**
     * Returns the cell type of a certain cell.
     * @param i Cell position in x direction.
     * @param j Cell position in y direction.
     * @param k Cell position in z direction.
     * @return The cell type (for more details @see GeometryValues).
     */
    inline uint32_t getCellType(int i, int j, int k) {
        return geometryValues.at(IDXFLAG(i,j,k));
    }

private:
    int imax, jmax, kmax;
    std::vector<uint32_t> geometryValues;
};

/**
 * Creates the .bingeo file for the passed scenario programmatically.
 * @param scenarioName The name of the scenario to create a geometry file for.
 * @param geometryFilename The filename of the .bingeo file to write.
 * @param imax Number of cells in x direction inside of the domain.
 * @param jmax Number of cells in y direction inside of the domain.
 * @param kmax Number of cells in z direction inside of the domain.
 */
void generateScenario(
        const std::string &scenarioName, const std::string &geometryFilename, int imax, int jmax, int kmax);

#endif //CFD3D_GEOMETRYCREATOR_HPP
