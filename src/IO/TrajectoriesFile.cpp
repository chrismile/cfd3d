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
#include <fstream>
#include <iomanip>
#include "IO/BinaryStream.hpp"
#include "TrajectoriesFile.hpp"

bool writeTrajectoriesToObjFile(const std::string &filename, const Trajectories &trajectories) {
    std::ofstream file(filename.c_str(), std::ofstream::binary);
    if (!file.is_open()) {
        std::cerr << "Error in writeLinesToObjFile: File \"" + filename + "\" not found." << std::endl;
        return false;
    }

    // We want five digits in output file
    file << std::setprecision(5);

    // Index of the next point
    size_t objPointIndex = 1;

    size_t trajectoryFileIndex = 0;
    for (size_t trajectoryIndex = 0; trajectoryIndex < trajectories.size(); trajectoryIndex++) {
        const Trajectory &trajectory = trajectories.at(trajectoryIndex);
        size_t trajectorySize = trajectory.positions.size();
        if (trajectorySize < 2) {
            continue;
        }

        for (size_t i = 0; i < trajectorySize; i++) {
            const glm::vec3 &v = trajectory.positions.at(i);
            file << "v " << std::setprecision(5) << v.x << " " << v.y << " " << v.z << "\n";
            file << "vt ";// << std::setprecision(5) << trajectory.attributes.at(0).at(i) << "\n";
            for (size_t attributeIndex = 0; attributeIndex < trajectory.attributes.size(); attributeIndex++) {
                file << std::setprecision(5) << trajectory.attributes.at(attributeIndex).at(i);
                if (attributeIndex != trajectory.attributes.size()-1) {
                    file << " ";
                }
            }
            file << "\n";
        }

        file << "g line" << trajectoryFileIndex << "\n";
        file << "l ";
        for (size_t i = 1; i < trajectorySize+1; i++) {
            file << objPointIndex << " ";
            objPointIndex++;
        }
        file << "\n\n";
        trajectoryFileIndex++;
    }

    file.close();
    return true;
}


const uint32_t LINE_FILE_FORMAT_VERSION = 1u;

bool writeTrajectoriesToBinLinesFile(const std::string &filename, const Trajectories &trajectories) {
    std::ofstream file(filename.c_str(), std::ofstream::binary);
    if (!file.is_open()) {
        std::cerr << "Error in writeLinesToBinLinesFile: File \"" + filename + "\" not found." << std::endl;
        return false;
    }

    uint32_t numTrajectories = trajectories.size();
    uint32_t numAttributes = trajectories.size() == 0 ? 0 : trajectories.at(0).attributes.size();

    BinaryWriteStream stream;
    stream.write((uint32_t)LINE_FILE_FORMAT_VERSION);
    stream.write((uint32_t)numTrajectories);
    stream.write((uint32_t)numAttributes);

    for (uint32_t trajectoryIndex = 0; trajectoryIndex < numTrajectories; trajectoryIndex++) {
        const Trajectory &currentTrajectory = trajectories.at(trajectoryIndex);
        size_t trajectoryNumPoints = currentTrajectory.positions.size();
        stream.write((void*)&currentTrajectory.positions.front(), sizeof(glm::vec3)*trajectoryNumPoints);

        assert(numAttributes == currentTrajectory.attributes.size());
        for (uint32_t attributeIndex = 0; attributeIndex < numAttributes; attributeIndex++) {
            const std::vector<float> &currentAttribute = currentTrajectory.attributes.at(attributeIndex);
            assert(trajectoryNumPoints == currentAttribute.size());
            stream.write((void*)&currentAttribute.front(), sizeof(float)*trajectoryNumPoints);
        }
    }

    file.write((const char*)stream.getBuffer(), stream.getSize());
    file.close();
    return true;
}
