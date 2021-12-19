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
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include "ScenarioFile.hpp"

std::map<std::string, std::string> loadVariablesFromDatFile(const std::string &filename) {
    std::map<std::string, std::string> variables;
    std::ifstream file(filename.c_str());
    if (!file.is_open()) {
        std::cerr << "Error in readScenarioConfigurationFromFile: File \"" << filename
                  << "\" does not exist." << std::endl;
        return variables;
    }

    std::string lineString;
    while (getline(file, lineString)) {
        while (!lineString.empty() && (lineString[lineString.size()-1] == '\r' || lineString[lineString.size()-1] == ' ')) {
            // Remove '\r' of Windows line ending
            lineString = lineString.substr(0, lineString.size() - 1);
        }

        // Skip empty lines or lines with a comment character.
        if (lineString.empty() || lineString.at(0) == '#') {
            continue;
        }

        // Split the line by spaces and tabulators.
        std::vector<std::string> line;
        boost::algorithm::split(line, lineString, boost::is_any_of("\t "), boost::token_compress_on);

        // Empty line with only whitespace characters.
        if (line.empty()) {
            continue;
        }
        assert(line.size() == 2);

        std::string key = line.at(0);
        std::string value = line.at(1);
        variables.insert(std::make_pair(key, value));
    }

    file.close();
    return variables;
}

std::string readStringVariable(const std::map<std::string, std::string> &variables, const std::string &name) {
    auto it = variables.find(name);
    if (it == variables.end()) {
        std::cerr << "Variable '" << name << "' not found in .dat file." << std::endl;
        exit(1);
    }
    return it->second;
}

std::string readStringVariableOptional(
        const std::map<std::string, std::string> &variables, const std::string &name,
        const std::string &fallback, bool &variableFound) {
    auto it = variables.find(name);
    if (it == variables.end()) {
        variableFound = false;
        return fallback;
    }
    return it->second;
}

int readIntVariable(const std::map<std::string, std::string> &variables, const std::string &name) {
    auto it = variables.find(name);
    if (it == variables.end()) {
        std::cerr << "Variable '" << name << "' not found in .dat file." << std::endl;
        exit(1);
    }
    return std::stoi(it->second);
}

int readIntVariableOptional(
        const std::map<std::string, std::string> &variables, const std::string &name,
        int fallback, bool &variableFound) {
    auto it = variables.find(name);
    if (it == variables.end()) {
        variableFound = false;
        return fallback;
    }
    return std::stoi(it->second);
}

Real readRealVariable(const std::map<std::string, std::string> &variables, const std::string &name) {
    auto it = variables.find(name);
    if (it == variables.end()) {
        std::cerr << "Variable '" << name << "' not found in .dat file." << std::endl;
        exit(1);
    }
    return stringToReal(it->second);
}

Real readRealVariableOptional(
        const std::map<std::string, std::string> &variables, const std::string &name,
        Real fallback, bool &variableFound) {
    auto it = variables.find(name);
    if (it == variables.end()) {
        variableFound = false;
        return fallback;
    }
    return stringToReal(it->second);
}

void readScenarioConfigurationFromFile(
        const std::string &scenarioFilename, std::string &scenarioName, std::string &geometryName,
        Real &tEnd, Real &dtWrite, Real &xLength, Real &yLength, Real &zLength,
        Real &xOrigin, Real &yOrigin, Real &zOrigin,
        Real &UI, Real &VI, Real &WI, Real &PI, Real &TI, Real &GX, Real &GY, Real &GZ,
        Real &Re, Real &Pr, Real &omg, Real &eps, int &itermax, Real &alpha, Real &beta, Real &dt, Real &tau,
        bool &useTemperature, Real &T_h, Real &T_c,
        int &imax, int &jmax, int &kmax, Real &dx, Real &dy, Real &dz) {
    std::map<std::string, std::string> variables = loadVariablesFromDatFile(scenarioFilename);

    scenarioName = readStringVariable(variables, "scenario");
    geometryName = readStringVariable(variables, "geometry");

    tEnd = readRealVariable(variables, "tEnd");
    dtWrite = readRealVariable(variables, "dtWrite");

    xLength = readRealVariable(variables, "xLength");
    yLength = readRealVariable(variables, "yLength");
    zLength = readRealVariable(variables, "zLength");
    xOrigin = readRealVariable(variables, "xOrigin");
    yOrigin = readRealVariable(variables, "yOrigin");
    zOrigin = readRealVariable(variables, "zOrigin");

    UI = readRealVariable(variables, "UI");
    VI = readRealVariable(variables, "VI");
    WI = readRealVariable(variables, "WI");
    PI = readRealVariable(variables, "PI");
    TI = readRealVariableOptional(variables, "TI", 0.0, useTemperature);
    GX = readRealVariable(variables, "GX");
    GY = readRealVariable(variables, "GY");
    GZ = readRealVariable(variables, "GZ");

    Re = readRealVariable(variables, "Re");
    Pr = readRealVariableOptional(variables, "Pr", 0.0, useTemperature);
    omg = readRealVariable(variables, "omg");
    eps = readRealVariable(variables, "eps");
    itermax = readIntVariable(variables, "itermax");
    alpha = readRealVariable(variables, "alpha");
    beta = readRealVariableOptional(variables, "beta", 0.0, useTemperature);
    dt = readRealVariable(variables, "dt");
    tau = readRealVariable(variables, "tau");

    T_c = readRealVariableOptional(variables, "T_c", 0.0, useTemperature);
    T_h = readRealVariableOptional(variables, "T_h", 0.0, useTemperature);

    imax = readIntVariable(variables, "imax");
    jmax = readIntVariable(variables, "jmax");
    kmax = readIntVariable(variables, "kmax");
    dx = xLength / Real(imax);
    dy = yLength / Real(jmax);
    dz = zLength / Real(kmax);
}
