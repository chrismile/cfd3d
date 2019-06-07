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
#include <chrono>
#include "CfdSolver/Init.hpp"
#include "CfdSolver/Flag.hpp"
#include "CfdSolver/CfdSolver.hpp"
#include "CfdSolver/Cpp/CfdSolverCpp.hpp"
#ifdef USE_CUDA
#include "CfdSolver/Cuda/CfdSolverCuda.hpp"
#endif
#ifdef USE_SYCL
#include "CfdSolver/Sycl/CfdSolverSycl.hpp"
#endif
#include "IO/IOUtils.hpp"
#include "IO/ArgumentParser.hpp"
#include "IO/ProgressBar.hpp"
#include "IO/ScenarioFile.hpp"
#include "IO/NetCdfWriter.hpp"
#include "IO/TrajectoriesFile.hpp"
#include "ParticleTracer/StreamlineTracer.hpp"
#include "ParticleTracer/StreaklineTracer.hpp"
#include "ParticleTracer/PathlineTracer.hpp"

const std::string outputDirectory = "output/";
const std::string scenarioDirectory = "../scenarios/";
const std::string geometryDirectory = "../geometry/";
const std::string lineDirectory = "lines/";

int main(int argc, char *argv[]) {
    CfdSolver *cfdSolver;
    ProgressBar progressBar;
    NetCdfWriter netCdfWriter;
    bool traceStreamlines = false, traceStreaklines = false, tracePathlines = false;
    std::vector<rvec3> particleSeedingLocations;
    bool dataIsUpToDate = true;
    bool shallWriteOutput = true;

    int imax, jmax, kmax, itermax;
    Real Re, Pr, UI, VI, WI, PI, TI, GX, GY, GZ, tEnd, dtWrite, xLength, yLength, zLength, xOrigin, yOrigin, zOrigin,
            dt, dx, dy, dz, alpha, omg, tau, eps, beta, T_h, T_c;
    bool useTemperature = true;
    std::string scenarioName, geometryName, scenarioFilename, geometryFilename, outputFilename, solverName;
    parseArguments(argc, argv, scenarioName, solverName);
    scenarioFilename = scenarioDirectory + scenarioName + ".dat";

    readScenarioConfigurationFromFile(
            scenarioFilename, scenarioName, geometryName,
            tEnd, dtWrite, xLength, yLength, zLength, xOrigin, yOrigin, zOrigin,
            UI, VI, WI, PI, TI, GX, GY, GZ,
            Re, Pr, omg, eps, itermax, alpha, beta, dt, tau, useTemperature,
            T_h, T_c, imax, jmax, kmax, dx, dy, dz);
    rvec3 gridOrigin = rvec3(xOrigin, yOrigin, zOrigin);
    rvec3 gridSize = rvec3(xLength, yLength, zLength);
    StreamlineTracer streamlineTracer;
    StreaklineTracer streaklineTracer(dtWrite*0.1); // Tracing frequency multiple of write time step.
    PathlineTracer pathlineTracer;

    if (!useTemperature){
        T_c = 0.0;
        T_h = 0.0;
        beta = 0.0;
        Pr = 0.0;
        TI = 0.0;
    }

    geometryFilename = geometryDirectory + geometryName;
    outputFilename = outputDirectory + scenarioName + ".nc";
    std::cout << "Scenario name: " << scenarioName << std::endl;
    std::cout << "Scenario file: " << scenarioFilename << std::endl;
    std::cout << "Geometry file: " << geometryFilename << std::endl;
    std::cout << "Output file: " << geometryFilename << std::endl;

    prepareOutputDirectory(outputDirectory, scenarioName);


    Real n = 0;
    Real t = 0;
    Real tWrite = 0;

    // Create all arrays for the simulation.
    Real *U = new Real[(imax+1)*(jmax+2)*(kmax+2)];
    Real *V = new Real[(imax+2)*(jmax+1)*(kmax+2)];
    Real *W = new Real[(imax+2)*(jmax+2)*(kmax+1)];
    Real *P = new Real[(imax+2)*(jmax+2)*(kmax+2)];
    Real *T = new Real[(imax+2)*(jmax+2)*(kmax+2)];
    FlagType *Flag = new FlagType[(imax+2)*(jmax+2)*(kmax+2)];

    if (geometryName == "none") {
        initFlagNoObstacles(scenarioName, imax, jmax, kmax, Flag);
    } else {
        initFlagFromGeometryFile(scenarioName, geometryFilename, imax, jmax, kmax, Flag);
    }
    initArrays(UI, VI, WI, PI, TI, imax, jmax, kmax, U, V, W, P, T, Flag);


    auto startTime = std::chrono::system_clock::now();

    if (shallWriteOutput) {
        netCdfWriter.openFile(outputFilename, imax, jmax, kmax, dx, dy, dz, xOrigin, yOrigin, zOrigin);
        netCdfWriter.writeTimestep(0, t, U, V, W, P, T, Flag);
    }


    if (solverName == "cpp") {
        cfdSolver = new CfdSolverCpp();
    }
#ifdef USE_CUDA
        else if (solverName == "cuda") {
        cfdSolver = new CfdSolverCuda();
    }
#endif
#ifdef USE_SYCL
        else if (solverName == "sycl") {
        cfdSolver = new CfdSolverSycl(imax, jmax, kmax, U, V, W, P, T, Flag);
    }
#endif
    else {
        std::cerr << "Fatal error: Unsupported solver name \"" << solverName << "\"." << std::endl;
    }
    cfdSolver->initialize(scenarioName, Re, Pr, omg, eps, itermax, alpha, beta, dt, tau, GX, GY, GZ, useTemperature,
            T_h, T_c, imax, jmax, kmax, dx, dy, dz, U, V, W, P, T, Flag);

    while (t < tEnd) {
        progressBar.printProgress(t, tEnd, 50);

        cfdSolver->setBoundaryValues();
        cfdSolver->setBoundaryValuesScenarioSpecific();
        dt = cfdSolver->calculateDt();
        if (useTemperature) {
            cfdSolver->calculateTemperature();
        }
        cfdSolver->calculateFgh();
        cfdSolver->calculateRs();

        cfdSolver->executeSorSolver();
        cfdSolver->calculateUvw();
        dataIsUpToDate = false;

        t += dt;
        tWrite += dt;
        n++;
        if (tWrite - dtWrite > 0) {
            if (shallWriteOutput) {
                if (!dataIsUpToDate) {
                    cfdSolver->getDataForOutput(U, V, W, P, T);
                    dataIsUpToDate = true;
                }
                netCdfWriter.writeTimestep(n, t, U, V, W, P, T, Flag);
            }
            progressBar.printOutput(n, t, 50);
            tWrite -= dtWrite;
        }
        if (traceStreaklines || tracePathlines) {
            if (!dataIsUpToDate) {
                cfdSolver->getDataForOutput(U, V, W, P, T);
                dataIsUpToDate = true;
            }
            if (traceStreaklines) {
                streaklineTracer.timeStep(t, dt, imax, jmax, kmax, dx, dy, dz, U, V, W, P, T);
            }
            if (tracePathlines) {
                pathlineTracer.timeStep(t, dt, imax, jmax, kmax, dx, dy, dz, U, V, W, P, T);
            }
        }
    }

    if (traceStreamlines) {
        if (!dataIsUpToDate) {
            cfdSolver->getDataForOutput(U, V, W, P, T);
            dataIsUpToDate = true;
        }
        Trajectories streamlines = streamlineTracer.trace(
                particleSeedingLocations, gridOrigin, gridSize, dt, imax, jmax, kmax, dx, dy, dz, U, V, W, P, T);
        writeTrajectoriesToObjFile(lineDirectory + scenarioName + "-streamlines.obj", streamlines);
    }

    if (traceStreaklines) {
        writeTrajectoriesToObjFile(lineDirectory + scenarioName + "-streaklines.obj", streaklineTracer.getTrajectories());
    }
    if (tracePathlines) {
        writeTrajectoriesToObjFile(lineDirectory + scenarioName + "-pathlines.obj", pathlineTracer.getTrajectories());
    }

    auto endTime = std::chrono::system_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    std::cout << "System time elapsed: " << (elapsedTime.count() * 1e-6) << "s" << std::endl;

    delete cfdSolver;
    delete[] U;
    delete[] V;
    delete[] W;
    delete[] P;
    delete[] T;
    delete[] Flag;

    return 0;
}