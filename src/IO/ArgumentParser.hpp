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

#ifndef CFD3D_ARGUMENTPARSER_HPP
#define CFD3D_ARGUMENTPARSER_HPP

#include <string>
#include "Defines.hpp"

class OutputFileWriter;

/**
 * Parses the command line arguments passed to the program. For more information on the format, please see README.md.
 * @param argc The number of arguments.
 * @param argv The arguments.
 * @param scenarioName The name of the scenario to use.
 * @param solverName The name of the solver to use.
 * @param outputFileWriterType The type of the output file writer to use.
 * @param shallWriteOutput Whether to write an output file at all.
 * @param numParticles The number of particles to seed when using a particle tracer.
 * @param traceStreamlines Whether to trace streamlines in the fluid flow.
 * @param traceStreaklines Whether to trace streaklines in the fluid flow.
 * @param tracePathlines Whether to trace pathlines in the fluid flow.
 * @param iproc The number of processes in x direction (MPI solver only).
 * @param jproc The number of processes in y direction (MPI solver only).
 * @param kproc The number of processes in z direction (MPI solver only).
 * @param blockSizeX The block size to use for 3D domains in x direction (CUDA solver only).
 * @param blockSizeY The block size to use for 3D domains in y direction (CUDA solver only).
 * @param blockSizeZ The block size to use for 3D domains in z direction (CUDA solver only).
 * @param blockSize1D The block size to use for 1D domains (CUDA solver only).
 */
void parseArguments(
        int argc, char *argv[], std::string &scenarioName, std::string &solverName,
        std::string &outputFileWriterType, bool &shallWriteOutput, LinearSystemSolverType &linearSystemSolverType,
        int &numParticles, bool &traceStreamlines, bool &traceStreaklines, bool &tracePathlines,
        int &iproc, int &jproc, int &kproc,
        int &blockSizeX, int &blockSizeY, int &blockSizeZ, int &blockSize1D);

#endif //CFD3D_ARGUMENTPARSER_HPP
