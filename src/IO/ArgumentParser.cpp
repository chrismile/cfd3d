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

#include <cstring>
#include "ArgumentParser.hpp"

void parseArguments(
        int argc, char *argv[], std::string &scenarioName, std::string &solverName, int &numParticles,
        bool &traceStreamlines,  bool &traceStreaklines,  bool &tracePathlines) {
    scenarioName = "driven_cavity";
    solverName = "cpp";
    numParticles = 1000;
    traceStreamlines = false;
    traceStreaklines = false;
    tracePathlines = false;

    // Go over command line arguments.
    for (int i = 1; i < argc; i += 2) {
        if (strcmp(argv[i], "--scenario") == 0 && i != argc-1) {
            scenarioName = argv[i+1];
        } else if (strcmp(argv[i], "--solver") == 0 && i != argc-1) {
            solverName = argv[i+1];
        } else if (strcmp(argv[i], "--numparticles") == 0 && i != argc-1) {
            numParticles = std::stoi(argv[i+1]);
        } else if (strcmp(argv[i], "--tracestreamlines") == 0 && i != argc-1) {
            traceStreamlines = strcmp(argv[i+1], "false") == 0 ? false : true;
        } else if (strcmp(argv[i], "--tracestreaklines") == 0 && i != argc-1) {
            traceStreaklines = strcmp(argv[i+1], "false") == 0 ? false : true;
        } else if (strcmp(argv[i], "--tracepathlines") == 0 && i != argc-1) {
            tracePathlines = strcmp(argv[i+1], "false") == 0 ? false : true;
        }
    }
}
