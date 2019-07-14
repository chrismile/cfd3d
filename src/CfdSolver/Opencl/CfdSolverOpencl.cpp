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

CfdSolverOpencl::CfdSolverOpencl(int platformId, int blockSizeX, int blockSizeY, int blockSizeZ, int blockSize1D) {
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
    queue = cl::CommandQueue(context, devices[0]);
#else
    queue = cl::CommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE);
#endif

    size_t maxWorkGroupSize;
    devices[0].getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &maxWorkGroupSize);

    // Set local work size.
    workGroupSize1D = cl::NDRange(blockSize1D);
    workGroupSize2D = cl::NDRange(blockSizeX, blockSizeY);
    workGroupSize3D = cl::NDRange(blockSizeX, blockSizeY, blockSizeZ);
    assert(blockSizeX * blockSizeY * blockSizeZ <= maxWorkGroupSize);
}

void CfdSolverOpencl::initialize(
        const std::string &scenarioName, LinearSystemSolverType linearSystemSolverType,
        Real Re, Real Pr, Real omg, Real eps, int itermax, Real alpha, Real beta, Real dt, Real tau,
        Real GX, Real GY, Real GZ, bool useTemperature, Real T_h, Real T_c,
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz,
        Real *U, Real *V, Real *W, Real *P, Real *T, uint32_t *Flag) {
    this->scenarioName = scenarioName;
    this->linearSystemSolverType = linearSystemSolverType;
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

    localMemoryReductionReal = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Real)*blockSize1D);
    localMemoryReductionUint = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(unsigned int)*blockSize1D);

    this->Up = new Real[(imax+1)*(jmax+2)*(kmax+2)];
    this->Vp = new Real[(imax+2)*(jmax+1)*(kmax+2)];
    this->Wp = new Real[(imax+2)*(jmax+2)*(kmax+1)];
    this->Pp = new Real[(imax+2)*(jmax+2)*(kmax+2)];
    this->P_tempp = new Real[(imax+2)*(jmax+2)*(kmax+2)];
    this->Tp = new Real[(imax+2)*(jmax+2)*(kmax+2)];
    this->T_tempp = new Real[(imax+2)*(jmax+2)*(kmax+2)];
    this->Fp = new Real[(imax+1)*(jmax+1)*(kmax+1)];
    this->Gp = new Real[(imax+1)*(jmax+1)*(kmax+1)];
    this->Hp = new Real[(imax+1)*(jmax+1)*(kmax+1)];
    this->RSp = new Real[(imax+1)*(jmax+1)*(kmax+1)];
    this->Flagp = new FlagType[(imax+2)*(jmax+2)*(kmax+2)];

}

CfdSolverOpencl::~CfdSolverOpencl() {
    delete[] Up;
    delete[] Vp;
    delete[] Wp;
    delete[] Pp;
    delete[] P_tempp;
    delete[] Tp;
    delete[] T_tempp;
    delete[] Fp;
    delete[] Gp;
    delete[] Hp;
    delete[] RSp;
    delete[] Flagp;
}

void CfdSolverOpencl::setBoundaryValues() {
    /*debugToCpu();
    setBoundaryValuesCpp(T_h, T_c, imax, jmax, kmax, Up, Vp, Wp, Tp, Flagp);
    debugToGpu();*/

    cl::EnqueueArgs eargsYZ(
            queue, cl::NullRange, ClInterface::get()->rangePadding2D(kmax, jmax, workGroupSize2D), workGroupSize2D);
    cl::make_kernel<Real, Real, int, int, int, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>
            setLeftRightBoundariesOpencl(setLeftRightBoundariesOpenclKernel);
    setLeftRightBoundariesOpencl(eargsYZ, T_h, T_c, imax, jmax, kmax, U, V, W, T, Flag);

    cl::EnqueueArgs eargsXZ(
            queue, cl::NullRange, ClInterface::get()->rangePadding2D(kmax, imax, workGroupSize2D), workGroupSize2D);
    cl::make_kernel<Real, Real, int, int, int, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>
            setDownUpBoundariesOpencl(setDownUpBoundariesOpenclKernel);
    setDownUpBoundariesOpencl(eargsXZ, T_h, T_c, imax, jmax, kmax, U, V, W, T, Flag);

    cl::EnqueueArgs eargsXY(
            queue, cl::NullRange, ClInterface::get()->rangePadding2D(jmax, imax, workGroupSize2D), workGroupSize2D);
    cl::make_kernel<Real, Real, int, int, int, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>
            setFrontBackBoundariesOpencl(setFrontBackBoundariesOpenclKernel);
    setFrontBackBoundariesOpencl(eargsXY, T_h, T_c, imax, jmax, kmax, U, V, W, T, Flag);


    cl::EnqueueArgs eargs3D(
            queue, cl::NullRange,
            ClInterface::get()->rangePadding3D(kmax, jmax, imax, workGroupSize3D), workGroupSize3D);
    cl::make_kernel<int, int, int, cl::Buffer, cl::Buffer>
            setInternalUBoundariesOpencl(setInternalUBoundariesOpenclKernel);
    setInternalUBoundariesOpencl(eargs3D, imax, jmax, kmax, U, Flag);
    cl::make_kernel<int, int, int, cl::Buffer, cl::Buffer>
            setInternalVBoundariesOpencl(setInternalVBoundariesOpenclKernel);
    setInternalVBoundariesOpencl(eargs3D, imax, jmax, kmax, V, Flag);
    cl::make_kernel<int, int, int, cl::Buffer, cl::Buffer>
            setInternalWBoundariesOpencl(setInternalWBoundariesOpenclKernel);
    setInternalWBoundariesOpencl(eargs3D, imax, jmax, kmax, W, Flag);
    cl::make_kernel<int, int, int, cl::Buffer, cl::Buffer>
            setInternalTBoundariesOpencl(setInternalTBoundariesOpenclKernel);
    setInternalTBoundariesOpencl(eargs3D, imax, jmax, kmax, T, Flag);
}

void CfdSolverOpencl::setBoundaryValuesScenarioSpecific() {
    /*debugToCpu();
    setBoundaryValuesScenarioSpecificCpp(scenarioName, imax, jmax, kmax, Up, Vp, Wp, Flagp);
    debugToGpu();*/

    if (scenarioName == "driven_cavity") {
        cl::EnqueueArgs eargsXZ(queue, cl::NullRange,
                ClInterface::get()->rangePadding2D(kmax, imax + 1, workGroupSize2D), workGroupSize2D);
        cl::make_kernel<int, int, int, cl::Buffer>
                setDrivenCavityBoundariesOpencl(setDrivenCavityBoundariesOpenclKernel);
        setDrivenCavityBoundariesOpencl(eargsXZ, imax, jmax, kmax, U);
    } else if (scenarioName == "flow_over_step") {
        cl::EnqueueArgs eargsYZ(queue, cl::NullRange,
                ClInterface::get()->rangePadding2D(kmax, jmax - (jmax / 2 + 1) + 1, workGroupSize2D), workGroupSize2D);
        cl::make_kernel<int, int, int, cl::Buffer, cl::Buffer, cl::Buffer>
                setFlowOverStepBoundariesOpencl(setFlowOverStepBoundariesOpenclKernel);
        setFlowOverStepBoundariesOpencl(eargsYZ, imax, jmax, kmax, U, V, W);
    } else if (scenarioName == "single_tower") {
        cl::EnqueueArgs eargsYZ(queue, cl::NullRange,
                ClInterface::get()->rangePadding2D(kmax, jmax, workGroupSize2D), workGroupSize2D);
        cl::make_kernel<int, int, int, cl::Buffer, cl::Buffer, cl::Buffer>
                setSingleTowerBoundariesOpencl(setSingleTowerBoundariesOpenclKernel);
        setSingleTowerBoundariesOpencl(eargsYZ, imax, jmax, kmax, U, V, W);
    } else if (scenarioName == "terrain_1" || scenarioName == "fuji_san" || scenarioName == "zugspitze") {
        cl::EnqueueArgs eargsYZ(queue, cl::NullRange,
                ClInterface::get()->rangePadding2D(kmax, jmax, workGroupSize2D), workGroupSize2D);
        cl::make_kernel<int, int, int, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>
                setMountainBoundariesOpencl(setMountainBoundariesOpenclKernel);
        setMountainBoundariesOpencl(eargsYZ, imax, jmax, kmax, U, V, W, Flag);
    }
}

Real CfdSolverOpencl::calculateDt() {
    /*debugToCpu();
    calculateDtCpp(Re, Pr, tau, dt, dx, dy, dz, imax, jmax, kmax, Up, Vp, Wp, useTemperature);
    debugToGpu();*/

    calculateDtOpencl(
            Re, Pr, tau, dt, dx, dy, dz, imax, jmax, kmax, blockSize1D,
            queue, workGroupSize1D, calculateMaximumKernel, U, V, W,
            openclReductionArrayU1, openclReductionArrayU2,
            openclReductionArrayV1, openclReductionArrayV2,
            openclReductionArrayW1, openclReductionArrayW2,
            localMemoryReductionReal, useTemperature);

    return dt;
}


void CfdSolverOpencl::calculateTemperature() {
    /*debugToCpu();
    Real *tempp = Tp;
    Tp = T_tempp;
    T_tempp = tempp;
    calculateTemperatureCpp(Re, Pr, alpha, dt, dx, dy, dz, imax, jmax, kmax, Up, Vp, Wp, Tp, T_tempp, Flagp);
    debugToGpu();*/

    cl::Buffer temp = T;
    T = T_temp;
    T_temp = temp;

    cl::EnqueueArgs eargs(
            queue, cl::NullRange,
            ClInterface::get()->rangePadding3D(kmax, jmax, imax, workGroupSize3D), workGroupSize3D);
    cl::make_kernel<Real, Real, Real, Real, Real, Real, Real, int, int, int,
        cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer> calculateTemperatureOpencl(
            calculateTemperatureOpenclKernel);
    calculateTemperatureOpencl(eargs, Re, Pr, alpha, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, T, T_temp, Flag);
}

void CfdSolverOpencl::calculateFgh() {
    /*debugToCpu();
    calculateFghCpp(Re, GX, GY, GZ, alpha, beta, dt, dx, dy, dz, imax, jmax, kmax, Up, Vp, Wp, Tp, Fp, Gp, Hp, Flagp);
    debugToGpu();*/

    cl::EnqueueArgs eargs3D(
            queue, cl::NullRange,
            ClInterface::get()->rangePadding3D(kmax, jmax, imax, workGroupSize3D), workGroupSize3D);
    cl::make_kernel<Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, int, int, int,
        cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>
            calculateFghOpencl(calculateFghOpenclKernel);
    calculateFghOpencl(
            eargs3D, Re, GX, GY, GZ, alpha, beta, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, T, F, G, H, Flag);

    cl::EnqueueArgs eargsYZ(queue, cl::NullRange,
            ClInterface::get()->rangePadding2D(kmax, jmax, workGroupSize2D), workGroupSize2D);
    cl::make_kernel<int, int, int, cl::Buffer, cl::Buffer> setFBoundariesOpencl(setFBoundariesOpenclKernel);
    setFBoundariesOpencl(eargsYZ, imax, jmax, kmax, U, F);

    cl::EnqueueArgs eargsXZ(queue, cl::NullRange,
            ClInterface::get()->rangePadding2D(kmax, imax, workGroupSize2D), workGroupSize2D);
    cl::make_kernel<int, int, int, cl::Buffer, cl::Buffer> setGBoundariesOpencl(setGBoundariesOpenclKernel);
    setGBoundariesOpencl(eargsXZ, imax, jmax, kmax, V, G);

    cl::EnqueueArgs eargsXY(queue, cl::NullRange,
            ClInterface::get()->rangePadding2D(jmax, imax, workGroupSize2D), workGroupSize2D);
    cl::make_kernel<int, int, int, cl::Buffer, cl::Buffer> setHBoundariesOpencl(setHBoundariesOpenclKernel);
    setHBoundariesOpencl(eargsXY, imax, jmax, kmax, W, H);
}

void CfdSolverOpencl::calculateRs() {
    /*debugToCpu();
    calculateRsCpp(dt, dx, dy, dz, imax, jmax, kmax, Fp, Gp, Hp, RSp);
    debugToGpu();*/

    cl::EnqueueArgs eargs(
            queue, cl::NullRange,
            ClInterface::get()->rangePadding3D(kmax, jmax, imax, workGroupSize3D), workGroupSize3D);
    cl::make_kernel<Real, Real, Real, Real, int, int, int, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>
        calculateRsOpencl(calculateRsOpenclKernel);
    calculateRsOpencl(eargs, dt, dx, dy, dz, imax, jmax, kmax, F, G, H, RS);
}


void CfdSolverOpencl::executeSorSolver() {
    /*debugToCpu();
    sorSolverCpp(omg, eps, itermax, linearSystemSolverType, dx, dy, dz, imax, jmax, kmax, Pp, P_tempp, RSp, Flagp);
    debugToGpu();*/

    sorSolverOpencl(
            omg, eps, itermax, linearSystemSolverType, dx, dy, dz, imax, jmax, kmax,
            blockSizeX, blockSizeY, blockSizeZ, blockSize1D,
            queue, workGroupSize1D, workGroupSize2D, workGroupSize3D,
            P, P_temp, RS, Flag,
            openclReductionArrayResidual1, openclReductionArrayResidual2,
            openclReductionArrayNumCells1, openclReductionArrayNumCells2,
            localMemoryReductionReal, localMemoryReductionUint,
            setXYPlanesPressureBoundariesOpenclKernel, setXZPlanesPressureBoundariesOpenclKernel,
            setYZPlanesPressureBoundariesOpenclKernel,
            setBoundaryConditionsPressureInDomainOpenclKernel, copyPressureOpenclKernel,
            reduceSumOpenclKernelReal, reduceSumOpenclKernelUint,
            sorSolverIterationOpenclKernel, sorSolverComputeResidualArrayOpenclKernel);
}

void CfdSolverOpencl::calculateUvw() {
    /*debugToCpu();
    calculateUvwCpp(dt, dx, dy, dz, imax, jmax, kmax, Up, Vp, Wp, Fp, Gp, Hp, Pp, Flagp);
    debugToGpu();*/


    cl::EnqueueArgs eargs(
            queue, cl::NullRange,
            ClInterface::get()->rangePadding3D(kmax, jmax, imax, workGroupSize3D), workGroupSize3D);
    cl::make_kernel<Real, Real, Real, Real, int, int, int, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer,
        cl::Buffer, cl::Buffer, cl::Buffer> calculateUvwOpencl(calculateUvwOpenclKernel);
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

void CfdSolverOpencl::debugToCpu()
{
    queue.enqueueReadBuffer(this->U, CL_FALSE, 0, sizeof(Real)*(imax+1)*(jmax+2)*(kmax+2), (void*)Up);
    queue.enqueueReadBuffer(this->V, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+1)*(kmax+2), (void*)Vp);
    queue.enqueueReadBuffer(this->W, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+1), (void*)Wp);
    queue.enqueueReadBuffer(this->P, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), (void*)Pp);
    queue.enqueueReadBuffer(this->P_temp, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), (void*)P_tempp);
    queue.enqueueReadBuffer(this->T, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), (void*)Tp);
    queue.enqueueReadBuffer(this->T_temp, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), (void*)T_tempp);
    queue.enqueueReadBuffer(this->F, CL_FALSE, 0, sizeof(Real)*(imax+1)*(jmax+1)*(kmax+1), (void*)Fp);
    queue.enqueueReadBuffer(this->G, CL_FALSE, 0, sizeof(Real)*(imax+1)*(jmax+1)*(kmax+1), (void*)Gp);
    queue.enqueueReadBuffer(this->H, CL_FALSE, 0, sizeof(Real)*(imax+1)*(jmax+1)*(kmax+1), (void*)Hp);
    queue.enqueueReadBuffer(this->RS, CL_FALSE, 0, sizeof(Real)*(imax+1)*(jmax+1)*(kmax+1), (void*)RSp);
    queue.enqueueReadBuffer(this->Flag, CL_FALSE, 0, sizeof(FlagType)*(imax+2)*(jmax+2)*(kmax+2), (void*)Flagp);
    queue.finish();
}

void CfdSolverOpencl::debugToGpu()
{
    queue.enqueueWriteBuffer(this->U, CL_FALSE, 0, sizeof(Real)*(imax+1)*(jmax+2)*(kmax+2), (void*)Up);
    queue.enqueueWriteBuffer(this->V, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+1)*(kmax+2), (void*)Vp);
    queue.enqueueWriteBuffer(this->W, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+1), (void*)Wp);
    queue.enqueueWriteBuffer(this->P, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), (void*)Pp);
    queue.enqueueWriteBuffer(this->P_temp, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), (void*)P_tempp);
    queue.enqueueWriteBuffer(this->T, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), (void*)Tp);
    queue.enqueueWriteBuffer(this->T_temp, CL_FALSE, 0, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2), (void*)T_tempp);
    queue.enqueueWriteBuffer(this->F, CL_FALSE, 0, sizeof(Real)*(imax+1)*(jmax+1)*(kmax+1), (void*)Fp);
    queue.enqueueWriteBuffer(this->G, CL_FALSE, 0, sizeof(Real)*(imax+1)*(jmax+1)*(kmax+1), (void*)Gp);
    queue.enqueueWriteBuffer(this->H, CL_FALSE, 0, sizeof(Real)*(imax+1)*(jmax+1)*(kmax+1), (void*)Hp);
    queue.enqueueWriteBuffer(this->RS, CL_FALSE, 0, sizeof(Real)*(imax+1)*(jmax+1)*(kmax+1), (void*)RSp);
    queue.finish();
}