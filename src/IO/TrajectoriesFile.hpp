//
// Created by christoph on 03.06.19.
//

#ifndef CFD3D_TRAJECTORIESFILE_HPP
#define CFD3D_TRAJECTORIESFILE_HPP

#include <string>
#include <vector>
#include "Defines.hpp"
#include "ParticleTracer/ParticleTracer.hpp"

/**
 * Saves a list of trajectories to an .obj file.
 * @param filename The file name of the .obj file.
 * @param trajectories The list of trajectories.
 * @return true, if the file could be opened for writing, otherwise false.
 */
bool writeTrajectoriesToObjFile(const std::string &filename, const Trajectories &trajectories);

/**
 * Saves a list of trajectories to a .binlines file.
 * @param filename The file name of the .obj file.
 * @param trajectories The list of trajectories.
 * @return true, if the file could be opened for writing, otherwise false.
 */
bool writeTrajectoriesToBinLinesFile(const std::string &filename, const Trajectories &trajectories);

#endif //CFD3D_TRAJECTORIESFILE_HPP
