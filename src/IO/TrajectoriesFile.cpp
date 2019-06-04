//
// Created by christoph on 03.06.19.
//

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
            file << "vt " << std::setprecision(5) << trajectory.attributes.at(0).at(i) << "\n";
            for (size_t attributeIndex = 0; attributeIndex < trajectory.attributes.size(); attributeIndex++) {
                file << std::setprecision(5) << trajectory.attributes.at(0).at(i);
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
