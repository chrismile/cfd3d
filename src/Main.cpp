//
// Created by christoph on 02.06.19.
//

#include <iostream>
#include <chrono>
#include "CfdSolver/Init.hpp"
#include "CfdSolver/Flag.hpp"
#include "CfdSolver/CfdSolver.hpp"
#include "CfdSolver/Cpp/CfdSolverCpp.hpp"
#include "IO/IOUtils.hpp"
#include "IO/ProgressBar.hpp"
#include "IO/ScenarioFile.hpp"
#include "IO/NetCdfWriter.hpp"
#include "IO/TrajectoriesFile.hpp"
#include "ParticleTracer/StreamlineTracer.hpp"
#include "ParticleTracer/StreaklineTracer.hpp"
#include "ParticleTracer/PathlineTracer.hpp"

const std::string outputDirectory = "output/";
const std::string scenarioDirectory = "scenarios/";
const std::string geometryDirectory = "geometry/";
const std::string lineDirectory = "lines/";

int main(int argc, char *argv[]) {
    CfdSolver *cfdSolver = new CfdSolverCpp();
    ProgressBar progressBar;
    NetCdfWriter netCdfWriter;
    bool traceStreamlines = false, traceStreaklines = false, tracePathlines = false;
    std::vector<rvec3> particleSeedingLocations;
    bool dataIsUpToDate = true;

    int imax, jmax, kmax, itermax;
    Real Re, Pr, UI, VI, WI, PI, TI, GX, GY, GZ, tEnd, dtWrite, xLength, yLength, zLength, xOrigin, yOrigin, zOrigin,
            dt, dx, dy, dz, alpha, omg, tau, eps, res, beta, T_h, T_c;
    bool useTemperature = true;
    std::string scenarioName, geometryName, scenarioFilename, geometryFilename, outputFilename;

    scenarioFilename = outputDirectory + "driven-cavity.dat";
    if (argc > 1) {
        scenarioFilename = scenarioDirectory + argv[1] + ".dat";
    }

    readScenarioConfigurationFromFile(scenarioFilename, scenarioName, geometryName,
            tEnd, dtWrite, xLength, yLength, zLength, xOrigin, yOrigin, zOrigin,
            Re, Pr, omg, eps, itermax, alpha, beta, dt, tau, GX, GY, GZ, useTemperature,
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

    netCdfWriter.openFile(outputFilename, imax, jmax, kmax, dx, dy, dz, xOrigin, yOrigin, zOrigin);
    netCdfWriter.writeTimestep(0, t, U, V, W, P, T, Flag);
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
            if (!dataIsUpToDate) {
                cfdSolver->getDataForOutput(U, V, W, P, T);
                dataIsUpToDate = true;
            }
            netCdfWriter.writeTimestep(n, t, U, V, W, P, T, Flag);
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
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "System time elapsed: " << (elapsedTime.count() * 1e-6) << std::endl << "s" << std::endl;

    delete cfdSolver;
    delete[] U;
    delete[] V;
    delete[] W;
    delete[] P;
    delete[] T;
    delete[] Flag;

    return 0;
}