//
// Created by christoph on 01.06.19.
//

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
        int imax, int jmax, int kmax, FlagType *&FlagPtr);

/**
 * For scenarios not using arbitrary geometries.
 * This function justs sets the whole interior domain to fluid cells and all boundary cells to no-slip cells.
 * @param scenarioName The name of the CFD scenario to simulate.
 * @param imax Number of cells in x direction inside of the domain.
 * @param jmax Number of cells in y direction inside of the domain.
 * @param kmax Number of cells in z direction inside of the domain.
 * @param FlagPtr A reference to a Flag array. This function uses new[] to allocate the memory of the Flag array.
 */
void initFlagNoObstacles(const std::string &scenarioName, int imax, int jmax, int kmax, FlagType *&FlagPtr);

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
