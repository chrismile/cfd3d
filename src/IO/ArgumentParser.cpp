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
#include <iostream>
#include "ArgumentParser.hpp"

void parseArguments(
        int argc, char *argv[], std::string &scenarioName, std::string &solverName,
        std::string &outputFileWriterType, bool &shallWriteOutput, LinearSystemSolverType &linearSystemSolverType,
        int &numParticles, bool &traceStreamlines, int &iproc, int &jproc, int &kproc,
        int &blockSizeX, int &blockSizeY, int &blockSizeZ, int &blockSize1D,
        int &gpuId, int &openclPlatformId) {
    // driven_cavity, natural_convection, rayleigh_benard_convection_8-2-1, flow_over_step, single_tower, terrain_1,
    // fuji_san, zugspitze, ...
    scenarioName = "driven_cavity";
    solverName = "cpp";
    shallWriteOutput = true;
    outputFileWriterType = "vtk";
    linearSystemSolverType = LINEAR_SOLVER_JACOBI;
    numParticles = 500;
    traceStreamlines = false;
    iproc = jproc = kproc = 1;
    blockSizeX = blockSizeY = 8;
    blockSizeZ = 4;
    blockSize1D = blockSizeX * blockSizeY * blockSizeZ;
    openclPlatformId = 0;

    // Go over command line arguments.
    for (int i = 1; i < argc; i += 2) {
        if (strcmp(argv[i], "--scenario") == 0 && i != argc - 1) {
            scenarioName = argv[i+1];
        } else if (strcmp(argv[i], "--solver") == 0 && i != argc - 1) {
            solverName = argv[i+1];
        } else if (strcmp(argv[i], "--outputformat") == 0 && i != argc - 1) {
            outputFileWriterType = argv[i+1];
        } else if (strcmp(argv[i], "--numparticles") == 0 && i != argc - 1) {
            numParticles = std::stoi(argv[i+1]);
        } else if (strcmp(argv[i], "--linsolver") == 0 && i != argc - 1) {
            if (strcmp(argv[i+1], "jacobi") == 0) {
                linearSystemSolverType = LINEAR_SOLVER_JACOBI;
            } else if (strcmp(argv[i+1], "sor") == 0) {
                linearSystemSolverType = LINEAR_SOLVER_SOR;
            } else if (strcmp(argv[i+1], "gauss-seidel") == 0) {
                linearSystemSolverType = LINEAR_SOLVER_GAUSS_SEIDEL;
            } else {
                std::cerr << "Specified invalid linear systems solver name." << std::endl;
                exit(1);
            }
            shallWriteOutput = strcmp(argv[i + 1], "false") != 0;
        } else if (strcmp(argv[i], "--output") == 0 && i != argc - 1) {
            shallWriteOutput = strcmp(argv[i + 1], "false") != 0;
        } else if (strcmp(argv[i], "--tracestreamlines") == 0 && i != argc - 1) {
            traceStreamlines = strcmp(argv[i + 1], "false") != 0;
        } else if (strcmp(argv[i], "--numproc") == 0 && i < argc - 3) {
            iproc = std::stoi(argv[i+1]);
            jproc = std::stoi(argv[i+2]);
            kproc = std::stoi(argv[i+3]);
            i += 2;
        } else if (strcmp(argv[i], "--blocksize") == 0 && i < argc - 3) {
            blockSizeX = std::stoi(argv[i+1]);
            blockSizeY = std::stoi(argv[i+2]);
            blockSizeZ = std::stoi(argv[i+3]);
            blockSize1D = blockSizeX*blockSizeY*blockSizeZ;
        } else if (strcmp(argv[i], "--gpu") == 0 && i < argc - 1) {
            gpuId = std::stoi(argv[i+1]);
        } else if (strcmp(argv[i], "--platformid") == 0 && i < argc - 1) {
            openclPlatformId = std::stoi(argv[i+1]);
        }
    }
}
