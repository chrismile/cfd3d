//
// Created by christoph on 02.06.19.
//

#ifndef CFD3D_IOUTILS_HPP
#define CFD3D_IOUTILS_HPP

#include <string>

/**
 * Makes sure the directory 'outputDirectory' exists.
 * @param outputDirectory The output directory name.
 * @param scenarioName The name of the scenario to simulate.
 */
void prepareOutputDirectory(const std::string &outputDirectory, const std::string &scenarioName);

#endif //CFD3D_IOUTILS_HPP
