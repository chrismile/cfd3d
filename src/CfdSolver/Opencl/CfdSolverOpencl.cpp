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
#include "UvwOpencl.hpp"
#include "SorSolverOpencl.hpp"
#include "CfdSolverOpencl.hpp"
#include "../../Defines.hpp"

#include "../Cpp/UvwCpp.hpp"
#include "../Cpp/SorSolverCpp.hpp"
#include "../Cpp/BoundaryValuesCpp.hpp"

const std::string openclSourceDirectory = "../src/CfdSolver/Opencl/";

CfdSolverOpencl::CfdSolverOpencl(int gpuId, int platformId, int blockSizeX, int blockSizeY, int blockSizeZ, int blockSize1D) {
    this->gpuId = gpuId;
    this->platformId = platformId;
    this->blockSizeX = blockSizeX;
    this->blockSizeY = blockSizeY;
    this->blockSizeZ = blockSizeZ;
    this->blockSize1D = blockSize1D;

    // Initialize the OpenCL context.
    ClInterface::get()->initialize(CLContextInfo(platformId));
    ClInterface::get()->printInfo();

    std::string prependHeader = "";
#ifdef REAL_FLOAT
    prependHeader += "#define REAL_FLOAT\n";
#else
    prependHeader += "#define REAL_DOUBLE\n";
#endif

    prependHeader += loadTextFile(openclSourceDirectory + "OpenclDefines.hpp");

    context = ClInterface::get()->getContext();
    devices = ClInterface::get()->getDevices();
    if (devices.empty()) {
        std::cerr << "Fatal error in CfdSolverOpencl::CfdSolverOpencl: No supported OpenCL devices were found."
                  << std::endl;
        exit(1);
    }
    if (gpuId >= int(devices.size())) {
        std::cerr << "Error in CfdSolverOpencl::CfdSolverOpencl: Invalid device ID specified. Setting device ID to 0."
                  << std::endl;
        gpuId = 0;
    }

    ClInterface::get()->setPrependHeader(prependHeader);
    computeProgramBoundaryValues = ClInterface::get()->loadProgramFromSourceFiles({
        openclSourceDirectory + "BoundaryValuesOpencl.cl"});
    computeProgramSor = ClInterface::get()->loadProgramFromSourceFiles({openclSourceDirectory + "SorSolverOpencl.cl"});
    computeProgramUvw = ClInterface::get()->loadProgramFromSourceFiles({openclSourceDirectory + "UvwOpencl.cl"});

    setLeftRightBoundariesOpenclKernel = cl::Kernel(computeProgramBoundaryValues, "setLeftRightBoundariesOpencl");
    setDownUpBoundariesOpenclKernel = cl::Kernel(computeProgramBoundaryValues, "setDownUpBoundariesOpencl");
    setFrontBackBoundariesOpenclKernel = cl::Kernel(computeProgramBoundaryValues, "setFrontBackBoundariesOpencl");
    setInternalUBoundariesOpenclKernel = cl::Kernel(computeProgramBoundaryValues, "setInternalUBoundariesOpencl");
    setInternalVBoundariesOpenclKernel = cl::Kernel(computeProgramBoundaryValues, "setInternalVBoundariesOpencl");
    setInternalWBoundariesOpenclKernel = cl::Kernel(computeProgramBoundaryValues, "setInternalWBoundariesOpencl");
    setInternalTBoundariesOpenclKernel = cl::Kernel(computeProgramBoundaryValues, "setInternalTBoundariesOpencl");
    setDrivenCavityBoundariesOpenclKernel = cl::Kernel(computeProgramBoundaryValues, "setDrivenCavityBoundariesOpencl");
    setFlowOverStepBoundariesOpenclKernel = cl::Kernel(computeProgramBoundaryValues, "setFlowOverStepBoundariesOpencl");
    setSingleTowerBoundariesOpenclKernel = cl::Kernel(computeProgramBoundaryValues, "setSingleTowerBoundariesOpencl");
    setMountainBoundariesOpenclKernel = cl::Kernel(computeProgramBoundaryValues, "setMountainBoundariesOpencl");

    setXYPlanesPressureBoundariesOpenclKernel = cl::Kernel(computeProgramSor, "setXYPlanesPressureBoundariesOpencl");
    setXZPlanesPressureBoundariesOpenclKernel = cl::Kernel(computeProgramSor, "setXZPlanesPressureBoundariesOpencl");
    setYZPlanesPressureBoundariesOpenclKernel = cl::Kernel(computeProgramSor, "setYZPlanesPressureBoundariesOpencl");
    setBoundaryConditionsPressureInDomainOpenclKernel = cl::Kernel(
            computeProgramSor, "setBoundaryConditionsPressureInDomainOpencl");
    copyPressureOpenclKernel = cl::Kernel(computeProgramSor, "copyPressureOpencl");
    reduceSumOpenclKernelReal = cl::Kernel(computeProgramSor, "reduceSumOpenclKernelReal");
    reduceSumOpenclKernelUint = cl::Kernel(computeProgramSor, "reduceSumOpenclKernelUint");
    sorSolverIterationOpenclKernel = cl::Kernel(computeProgramSor, "sorSolverIterationOpencl");
    sorSolverComputeResidualArrayOpenclKernel = cl::Kernel(computeProgramSor, "sorSolverComputeResidualArrayOpencl");

    calculateFghOpenclKernel = cl::Kernel(computeProgramUvw, "calculateFghOpencl");
    setFBoundariesOpenclKernel = cl::Kernel(computeProgramUvw, "setFBoundariesOpencl");
    setGBoundariesOpenclKernel = cl::Kernel(computeProgramUvw, "setGBoundariesOpencl");
    setHBoundariesOpenclKernel = cl::Kernel(computeProgramUvw, "setHBoundariesOpencl");
    calculateRsOpenclKernel = cl::Kernel(computeProgramUvw, "calculateRsOpencl");
    calculateMaximumKernel = cl::Kernel(computeProgramUvw, "calculateMaximum");
    calculateUvwOpenclKernel = cl::Kernel(computeProgramUvw, "calculateUvwOpencl");
    calculateTemperatureOpenclKernel = cl::Kernel(computeProgramUvw, "calculateTemperatureOpencl");

#ifndef _PROFILING_CL_
    queue = cl::CommandQueue(context, devices[gpuId]);
#else
    queue = cl::CommandQueue(context, devices[gpuId], CL_QUEUE_PROFILING_ENABLE);
#endif

    size_t maxWorkGroupSize;
    devices[gpuId].getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &maxWorkGroupSize);

    // Set local work size.
    workGroupSize1D = cl::NDRange(blockSize1D);
    workGroupSize2D = cl::NDRange(blockSizeX, blockSizeY);
    workGroupSize3D = cl::NDRange(blockSizeX, blockSizeY, blockSizeZ);
    assert(size_t(blockSizeX) * size_t(blockSizeY) * size_t(blockSizeZ) <= maxWorkGroupSize);
}

void CfdSolverOpencl::initialize(
        const std::string &scenarioName, LinearSystemSolverType linearSystemSolverType, bool shallWriteOutput,
        Real Re, Real Pr, Real omg, Real eps, int itermax, Real alpha, Real beta, Real dt, Real tau,
        Real GX, Real GY, Real GZ, bool useTemperature, Real T_h, Real T_c,
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz,
        Real *U, Real *V, Real *W, Real *P, Real *T, uint32_t *Flag) {
    this->scenarioName = scenarioName;
    this->linearSystemSolverType = linearSystemSolverType;
    this->shallWriteOutput = shallWriteOutput;
    this->Re = Re;
    this->Pr = Pr;
    this->omg = omg;
    this->eps = eps;
    this->itermax = itermax;
    this->alpha = alpha;
    this->beta = beta;
    this->dt = dt;
    this->tau = tau;
    this->GX = GX;
    this->GY = GY;
    this->GZ = GZ;
    this->useTemperature = useTemperature;
    this->T_h = T_h;
    this->T_c = T_c;
    this->imax = imax;
    this->jmax = jmax;
    this->kmax = kmax;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;

    // Create all arrays for the simulation.
    this->U = cl::Buffer(
            context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
            (imax+1)*(jmax+2)*(kmax+2)*sizeof(Real), (void*)U);
    this->V = cl::Buffer(
            context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
            (imax+2)*(jmax+1)*(kmax+2)*sizeof(Real), (void*)V);
    this->W = cl::Buffer(
            context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
            (imax+2)*(jmax+2)*(kmax+1)*sizeof(Real), (void*)W);
    this->P = cl::Buffer(
            context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
            (imax+2)*(jmax+2)*(kmax+2)*sizeof(Real), (void*)P);
    this->T = cl::Buffer(
            context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
            (imax+2)*(jmax+2)*(kmax+2)*sizeof(Real), (void*)T);
    this->Flag = cl::Buffer(
            context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            (imax+2)*(jmax+2)*(kmax+2)*sizeof(unsigned int), (void*)Flag);
    this->P_temp = cl::Buffer(context, CL_MEM_READ_WRITE, (imax+2)*(jmax+2)*(kmax+2)*sizeof(Real));
    this->T_temp = cl::Buffer(context, CL_MEM_READ_WRITE, (imax+2)*(jmax+2)*(kmax+2)*sizeof(Real));
    this->F = cl::Buffer(context, CL_MEM_READ_WRITE, (imax+1)*(jmax+1)*(kmax+1)*sizeof(Real));
    this->G = cl::Buffer(context, CL_MEM_READ_WRITE, (imax+1)*(jmax+1)*(kmax+1)*sizeof(Real));
    this->H = cl::Buffer(context, CL_MEM_READ_WRITE, (imax+1)*(jmax+1)*(kmax+1)*sizeof(Real));
    this->RS = cl::Buffer(context, CL_MEM_READ_WRITE, (imax+1)*(jmax+1)*(kmax+1)*sizeof(Real));

    int openclReductionArrayUSize = iceil((imax+1)*(jmax+2)*(kmax+2), blockSize1D*2);
    int openclReductionArrayVSize = iceil((imax+1)*(jmax+2)*(kmax+2), blockSize1D*2);
    int openclReductionArrayWSize = iceil((imax+1)*(jmax+2)*(kmax+2), blockSize1D*2);
    int openclReductionArrayResidualSize1 = iceil(imax*jmax*kmax, blockSize1D*2)*blockSize1D*2;
    int openclReductionArrayResidualSize2 = iceil(imax*jmax*kmax, blockSize1D*2);
    this->openclReductionArrayU1 = cl::Buffer(context, CL_MEM_READ_WRITE, openclReductionArrayUSize*sizeof(Real));
    this->openclReductionArrayU2 = cl::Buffer(context, CL_MEM_READ_WRITE, openclReductionArrayUSize*sizeof(Real));
    this->openclReductionArrayV1 = cl::Buffer(context, CL_MEM_READ_WRITE, openclReductionArrayVSize*sizeof(Real));
    this->openclReductionArrayV2 = cl::Buffer(context, CL_MEM_READ_WRITE, openclReductionArrayVSize*sizeof(Real));
    this->openclReductionArrayW1 = cl::Buffer(context, CL_MEM_READ_WRITE, openclReductionArrayWSize*sizeof(Real));
    this->openclReductionArrayW2 = cl::Buffer(context, CL_MEM_READ_WRITE, openclReductionArrayWSize*sizeof(Real));
    this->openclReductionArrayResidual1 = cl::Buffer(
            context, CL_MEM_READ_WRITE, openclReductionArrayResidualSize1*sizeof(Real));
    this->openclReductionArrayResidual2 = cl::Buffer(
            context, CL_MEM_READ_WRITE, openclReductionArrayResidualSize2*sizeof(Real));
    this->openclReductionArrayNumCells1 = cl::Buffer(
            context, CL_MEM_READ_WRITE, openclReductionArrayResidualSize1*sizeof(unsigned int));
    this->openclReductionArrayNumCells2 = cl::Buffer(
            context, CL_MEM_READ_WRITE, openclReductionArrayResidualSize2*sizeof(unsigned int));
}

CfdSolverOpencl::~CfdSolverOpencl() {
}

void CfdSolverOpencl::setBoundaryValues() {
    cl::EnqueueArgs eargsYZ(
            queue, cl::NullRange, ClInterface::get()->rangePadding2D(kmax, jmax, workGroupSize2D), workGroupSize2D);
    auto setLeftRightBoundariesOpencl =
            cl::KernelFunctor<Real, Real, int, int, int, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>(
                    setLeftRightBoundariesOpenclKernel);
    setLeftRightBoundariesOpencl(eargsYZ, T_h, T_c, imax, jmax, kmax, U, V, W, T, Flag);

    cl::EnqueueArgs eargsXZ(
            queue, cl::NullRange, ClInterface::get()->rangePadding2D(kmax, imax, workGroupSize2D), workGroupSize2D);
    auto setDownUpBoundariesOpencl =
            cl::KernelFunctor<Real, Real, int, int, int, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>(
                    setDownUpBoundariesOpenclKernel);
    setDownUpBoundariesOpencl(eargsXZ, T_h, T_c, imax, jmax, kmax, U, V, W, T, Flag);

    cl::EnqueueArgs eargsXY(
            queue, cl::NullRange, ClInterface::get()->rangePadding2D(jmax, imax, workGroupSize2D), workGroupSize2D);
    auto setFrontBackBoundariesOpencl =
            cl::KernelFunctor<Real, Real, int, int, int, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>(
                    setFrontBackBoundariesOpenclKernel);
    setFrontBackBoundariesOpencl(eargsXY, T_h, T_c, imax, jmax, kmax, U, V, W, T, Flag);


    cl::EnqueueArgs eargs3D(
            queue, cl::NullRange,
            ClInterface::get()->rangePadding3D(kmax, jmax, imax, workGroupSize3D), workGroupSize3D);
    auto setInternalUBoundariesOpencl = cl::KernelFunctor<int, int, int, cl::Buffer, cl::Buffer>(
            setInternalUBoundariesOpenclKernel);
    setInternalUBoundariesOpencl(eargs3D, imax, jmax, kmax, U, Flag);
    auto setInternalVBoundariesOpencl = cl::KernelFunctor<int, int, int, cl::Buffer, cl::Buffer>(
            setInternalVBoundariesOpenclKernel);
    setInternalVBoundariesOpencl(eargs3D, imax, jmax, kmax, V, Flag);
    auto setInternalWBoundariesOpencl = cl::KernelFunctor<int, int, int, cl::Buffer, cl::Buffer>(
            setInternalWBoundariesOpenclKernel);
    setInternalWBoundariesOpencl(eargs3D, imax, jmax, kmax, W, Flag);
    auto setInternalTBoundariesOpencl = cl::KernelFunctor<int, int, int, cl::Buffer, cl::Buffer>(
            setInternalTBoundariesOpenclKernel);
    setInternalTBoundariesOpencl(eargs3D, imax, jmax, kmax, T, Flag);
}

void CfdSolverOpencl::setBoundaryValuesScenarioSpecific() {
    if (scenarioName == "driven_cavity") {
        cl::EnqueueArgs eargsXZ(queue, cl::NullRange,
                ClInterface::get()->rangePadding2D(kmax, imax + 1, workGroupSize2D), workGroupSize2D);
        auto setDrivenCavityBoundariesOpencl = cl::KernelFunctor<int, int, int, cl::Buffer>(
                setDrivenCavityBoundariesOpenclKernel);
        setDrivenCavityBoundariesOpencl(eargsXZ, imax, jmax, kmax, U);
    } else if (scenarioName == "flow_over_step") {
        cl::EnqueueArgs eargsYZ(queue, cl::NullRange,
                ClInterface::get()->rangePadding2D(kmax, jmax - (jmax / 2 + 1) + 1, workGroupSize2D), workGroupSize2D);
        auto setFlowOverStepBoundariesOpencl = cl::KernelFunctor<int, int, int, cl::Buffer, cl::Buffer, cl::Buffer>(
                setFlowOverStepBoundariesOpenclKernel);
        setFlowOverStepBoundariesOpencl(eargsYZ, imax, jmax, kmax, U, V, W);
    } else if (scenarioName == "single_tower") {
        cl::EnqueueArgs eargsYZ(queue, cl::NullRange,
                ClInterface::get()->rangePadding2D(kmax, jmax, workGroupSize2D), workGroupSize2D);
        auto setSingleTowerBoundariesOpencl = cl::KernelFunctor<int, int, int, cl::Buffer, cl::Buffer, cl::Buffer>(
                setSingleTowerBoundariesOpenclKernel);
        setSingleTowerBoundariesOpencl(eargsYZ, imax, jmax, kmax, U, V, W);
    } else if (scenarioName == "terrain_1" || scenarioName == "fuji_san" || scenarioName == "zugspitze") {
        cl::EnqueueArgs eargsYZ(queue, cl::NullRange,
                ClInterface::get()->rangePadding2D(kmax, jmax, workGroupSize2D), workGroupSize2D);
        auto setMountainBoundariesOpencl =
                cl::KernelFunctor<int, int, int, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>(
                        setMountainBoundariesOpenclKernel);
        setMountainBoundariesOpencl(eargsYZ, imax, jmax, kmax, U, V, W, Flag);
    }
}

Real CfdSolverOpencl::calculateDt() {
    calculateDtOpencl(
            Re, Pr, tau, dt, dx, dy, dz, imax, jmax, kmax, blockSize1D,
            queue, workGroupSize1D, calculateMaximumKernel, U, V, W,
            openclReductionArrayU1, openclReductionArrayU2,
            openclReductionArrayV1, openclReductionArrayV2,
            openclReductionArrayW1, openclReductionArrayW2,
            useTemperature);

    return dt;
}


void CfdSolverOpencl::calculateTemperature() {
    cl::Buffer temp = T;
    T = T_temp;
    T_temp = temp;

    cl::EnqueueArgs eargs(
            queue, cl::NullRange,
            ClInterface::get()->rangePadding3D(kmax, jmax, imax, workGroupSize3D), workGroupSize3D);
    auto calculateTemperatureOpencl = cl::KernelFunctor<Real, Real, Real, Real, Real, Real, Real, int, int, int,
            cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>(
            calculateTemperatureOpenclKernel);
    calculateTemperatureOpencl(eargs, Re, Pr, alpha, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, T, T_temp, Flag);
}

void CfdSolverOpencl::calculateFgh() {
    cl::EnqueueArgs eargs3D(
            queue, cl::NullRange,
            ClInterface::get()->rangePadding3D(kmax, jmax, imax, workGroupSize3D), workGroupSize3D);
    auto calculateFghOpencl = cl::KernelFunctor<Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, int, int,
            int, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>(
            calculateFghOpenclKernel);
    calculateFghOpencl(
            eargs3D, Re, GX, GY, GZ, alpha, beta, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, T, F, G, H, Flag);

    cl::EnqueueArgs eargsYZ(queue, cl::NullRange,
            ClInterface::get()->rangePadding2D(kmax, jmax, workGroupSize2D), workGroupSize2D);
    auto setFBoundariesOpencl = cl::KernelFunctor<int, int, int, cl::Buffer, cl::Buffer>(setFBoundariesOpenclKernel);
    setFBoundariesOpencl(eargsYZ, imax, jmax, kmax, U, F);

    cl::EnqueueArgs eargsXZ(queue, cl::NullRange,
            ClInterface::get()->rangePadding2D(kmax, imax, workGroupSize2D), workGroupSize2D);
    auto setGBoundariesOpencl = cl::KernelFunctor<int, int, int, cl::Buffer, cl::Buffer>(setGBoundariesOpenclKernel);
    setGBoundariesOpencl(eargsXZ, imax, jmax, kmax, V, G);

    cl::EnqueueArgs eargsXY(queue, cl::NullRange,
            ClInterface::get()->rangePadding2D(jmax, imax, workGroupSize2D), workGroupSize2D);
    auto setHBoundariesOpencl = cl::KernelFunctor<int, int, int, cl::Buffer, cl::Buffer>(setHBoundariesOpenclKernel);
    setHBoundariesOpencl(eargsXY, imax, jmax, kmax, W, H);
}

void CfdSolverOpencl::calculateRs() {
    cl::EnqueueArgs eargs(
            queue, cl::NullRange,
            ClInterface::get()->rangePadding3D(kmax, jmax, imax, workGroupSize3D), workGroupSize3D);
    auto calculateRsOpencl =
            cl::KernelFunctor<Real, Real, Real, Real, int, int, int, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>(
                    calculateRsOpenclKernel);
    calculateRsOpencl(eargs, dt, dx, dy, dz, imax, jmax, kmax, F, G, H, RS);
}


void CfdSolverOpencl::executeSorSolver() {
    sorSolverOpencl(
            omg, eps, itermax, linearSystemSolverType, shallWriteOutput, dx, dy, dz, imax, jmax, kmax,
            blockSizeX, blockSizeY, blockSizeZ, blockSize1D,
            queue, workGroupSize1D, workGroupSize2D, workGroupSize3D,
            P, P_temp, RS, Flag,
            openclReductionArrayResidual1, openclReductionArrayResidual2,
            openclReductionArrayNumCells1, openclReductionArrayNumCells2,
            setXYPlanesPressureBoundariesOpenclKernel, setXZPlanesPressureBoundariesOpenclKernel,
            setYZPlanesPressureBoundariesOpenclKernel,
            setBoundaryConditionsPressureInDomainOpenclKernel, copyPressureOpenclKernel,
            reduceSumOpenclKernelReal, reduceSumOpenclKernelUint,
            sorSolverIterationOpenclKernel, sorSolverComputeResidualArrayOpenclKernel);
}

void CfdSolverOpencl::calculateUvw() {
    cl::EnqueueArgs eargs(
            queue, cl::NullRange,
            ClInterface::get()->rangePadding3D(kmax, jmax, imax, workGroupSize3D), workGroupSize3D);
    auto calculateUvwOpencl = cl::KernelFunctor<Real, Real, Real, Real, int, int, int, cl::Buffer, cl::Buffer,
            cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>(calculateUvwOpenclKernel);
    calculateUvwOpencl(eargs, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, F, G, H, P, Flag);
}

void CfdSolverOpencl::getDataForOutput(Real *U, Real *V, Real *W, Real *P, Real *T) {
    // Copy the content of U, V, W, P, T in the internal representation to the specified output arrays.
    queue.enqueueReadBuffer(this->U, CL_FALSE, 0, sizeof(Real)*(imax+1)*(jmax+2)*(kmax+2), (void*)U);
    queue.enqueueReadBuffer(this->V, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+1)*(kmax+2), (void*)V);
    queue.enqueueReadBuffer(this->W, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+1), (void*)W);
    queue.enqueueReadBuffer(this->P, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), (void*)P);
    queue.enqueueReadBuffer(this->T, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), (void*)T);
    queue.finish();
}
