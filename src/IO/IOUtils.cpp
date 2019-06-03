//
// Created by christoph on 02.06.19.
//

#include <iostream>
#include <cstdlib>
#include <boost/filesystem.hpp>
#include "IOUtils.hpp"

void prepareOutputDirectory(const std::string &outputDirectory, const std::string &scenarioName) {
    // Create the output directory.
    if (!boost::filesystem::exists(outputDirectory) && boost::filesystem::is_directory(outputDirectory)) {
        if (!boost::filesystem::create_directory(outputDirectory)) {
            std::cerr << "Output directory could not be created" << std::endl;
            exit(1);
        }
    }
}