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
#include <cmath>
#include "SorSolverOpencl.hpp"

/*
 * Overwrites the contents of the passed reduction arrays.
 */
Real reduceSumOpenclReal(
        cl::CommandQueue &queue, cl::NDRange workGroupSize1D, cl::Kernel &calculateSumKernel,
        cl::Buffer &input, unsigned int numValues, cl::Buffer &openclReductionHelperArray, int blockSize1D) {
    Real sumValue = Real(0);

    auto calculateSum = cl::KernelFunctor<cl::Buffer, cl::Buffer, cl::LocalSpaceArg, int>(calculateSumKernel);

    cl::Buffer *reductionInput = &input;
    cl::Buffer *reductionOutput = &openclReductionHelperArray;

    int numberOfBlocks = numValues;
    int inputSize;
    bool finished = false;

    int iteration = 0;
    while (!finished) {
        inputSize = numberOfBlocks;
        numberOfBlocks = iceil(numberOfBlocks, blockSize1D*2);

        if (inputSize != 1) {
            int localMemorySize = blockSize1D * sizeof(Real);
            cl::EnqueueArgs eargs(queue, cl::NullRange, cl::NDRange(numberOfBlocks*blockSize1D*2), workGroupSize1D);
            calculateSum(eargs, *reductionInput, *reductionOutput, cl::Local(localMemorySize), inputSize);
            if (iteration % 2 == 0) {
                reductionInput = &openclReductionHelperArray;
                reductionOutput = &input;
            } else {
                reductionInput = &input;
                reductionOutput = &openclReductionHelperArray;
            }
        }

        if (numberOfBlocks == 1) {
            finished = true;
        }
        iteration++;
    }

    queue.enqueueReadBuffer(*reductionInput, CL_FALSE, 0, sizeof(Real), (void*)&sumValue);
    queue.finish();

    return sumValue;
}

/*
 * Overwrites the contents of the passed reduction arrays.
 */
unsigned int reduceSumOpenclUint(
        cl::CommandQueue &queue, cl::NDRange workGroupSize1D, cl::Kernel &calculateSumKernel,
        cl::Buffer &input, unsigned int numValues, cl::Buffer &openclReductionHelperArray, int blockSize1D) {
    unsigned int sumValue = 0u;

    auto calculateSum = cl::KernelFunctor<cl::Buffer, cl::Buffer, cl::LocalSpaceArg, int>(calculateSumKernel);

    cl::Buffer *reductionInput = &input;
    cl::Buffer *reductionOutput = &openclReductionHelperArray;

    int numberOfBlocks = numValues;
    int inputSize;
    bool finished = false;

    int iteration = 0;
    while (!finished) {
        inputSize = numberOfBlocks;
        numberOfBlocks = iceil(numberOfBlocks, blockSize1D*2);

        if (inputSize != 1) {
            int localMemorySize = blockSize1D * sizeof(Real);
            cl::EnqueueArgs eargs(queue, cl::NullRange, cl::NDRange(numberOfBlocks*blockSize1D*2), workGroupSize1D);
            calculateSum(eargs, *reductionInput, *reductionOutput, cl::Local(localMemorySize), inputSize);
            if (iteration % 2 == 0) {
                reductionInput = &openclReductionHelperArray;
                reductionOutput = &input;
            } else {
                reductionInput = &input;
                reductionOutput = &openclReductionHelperArray;
            }
        }

        if (numberOfBlocks == 1) {
            finished = true;
        }
        iteration++;
    }

    queue.enqueueReadBuffer(*reductionInput, CL_FALSE, 0, sizeof(unsigned int), (void*)&sumValue);
    queue.finish();

    return sumValue;
}



void sorSolverOpencl(
        Real omg, Real eps, int itermax, LinearSystemSolverType linearSystemSolverType, bool shallWriteOutput,
        Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        int blockSizeX, int blockSizeY, int blockSizeZ, int blockSize1D,
        cl::CommandQueue &queue, cl::NDRange workGroupSize1D, cl::NDRange workGroupSize2D, cl::NDRange workGroupSize3D,
        cl::Buffer &P, cl::Buffer &P_temp, cl::Buffer &RS, cl::Buffer &Flag,
        cl::Buffer &openclReductionArrayResidual1, cl::Buffer &openclReductionArrayResidual2,
        cl::Buffer &openclReductionArrayNumCells1, cl::Buffer &openclReductionArrayNumCells2,
        cl::Kernel &setXYPlanesPressureBoundariesOpenclKernel, cl::Kernel &setXZPlanesPressureBoundariesOpenclKernel,
        cl::Kernel &setYZPlanesPressureBoundariesOpenclKernel,
        cl::Kernel &setBoundaryConditionsPressureInDomainOpenclKernel, cl::Kernel &copyPressureOpenclKernel,
        cl::Kernel &reduceSumOpenclKernelReal, cl::Kernel &reduceSumOpenclKernelUint,
        cl::Kernel &sorSolverIterationOpenclKernel, cl::Kernel &sorSolverComputeResidualArrayOpenclKernel) {
    if (linearSystemSolverType == LINEAR_SOLVER_SOR || linearSystemSolverType == LINEAR_SOLVER_SOR_PARALLEL) {
        // Successive over-relaxation based on Gauss-Seidl. A factor of 1.5 proved to give the best results here.
        omg = 1.5;
    } else {
        // A method named JOR (Jacobi over-relaxation) with omega != 1 exists, but doesn't converge for this problem.
        omg = 1.0;
    }


    cl::EnqueueArgs eargsXY(queue, cl::NullRange,
            ClInterface::get()->rangePadding2D(jmax, imax, workGroupSize2D), workGroupSize2D);
    auto setXYPlanesPressureBoundariesOpencl = cl::KernelFunctor<int, int, int, cl::Buffer>(
            setXYPlanesPressureBoundariesOpenclKernel);

    cl::EnqueueArgs eargsXZ(queue, cl::NullRange,
            ClInterface::get()->rangePadding2D(kmax, imax, workGroupSize2D), workGroupSize2D);
    auto setXZPlanesPressureBoundariesOpencl = cl::KernelFunctor<int, int, int, cl::Buffer>(
            setXZPlanesPressureBoundariesOpenclKernel);

    cl::EnqueueArgs eargsYZ(queue, cl::NullRange,
            ClInterface::get()->rangePadding2D(kmax, jmax, workGroupSize2D), workGroupSize2D);
    auto setYZPlanesPressureBoundariesOpencl = cl::KernelFunctor<int, int, int, cl::Buffer>(
            setYZPlanesPressureBoundariesOpenclKernel);

    cl::EnqueueArgs eargs3D(
            queue, cl::NullRange,
            ClInterface::get()->rangePadding3D(kmax, jmax, imax, workGroupSize3D), workGroupSize3D);
    auto setBoundaryConditionsPressureInDomainOpencl = cl::KernelFunctor<int, int, int, cl::Buffer, cl::Buffer>(
            setBoundaryConditionsPressureInDomainOpenclKernel);

    cl::EnqueueArgs eargsWholeDomain3D(
            queue, cl::NullRange,
            ClInterface::get()->rangePadding3D(kmax+2, jmax+2, imax+2, workGroupSize3D), workGroupSize3D);
    auto copyPressureOpencl = cl::KernelFunctor<int, int, int, cl::Buffer, cl::Buffer>(copyPressureOpenclKernel);

    auto sorSolverIterationOpencl = cl::KernelFunctor<
            Real, Real, Real, Real, Real, int, int, int, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>(
            sorSolverIterationOpenclKernel);

    auto sorSolverComputeResidualArrayOpencl = cl::KernelFunctor<
            Real, Real, Real, int, int, int, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer>(
            sorSolverComputeResidualArrayOpenclKernel);


    const Real coeff = Real(omg / (2.0 * (1.0 / (dx*dx) + 1.0 / (dy*dy) + 1.0 / (dz*dz))));
    Real residual = 1e9;
    int it = 0;
    while (it < itermax && residual > eps) {
        setXYPlanesPressureBoundariesOpencl(eargsXY, imax, jmax, kmax, P);
        setXZPlanesPressureBoundariesOpencl(eargsXZ, imax, jmax, kmax, P);
        setYZPlanesPressureBoundariesOpencl(eargsYZ, imax, jmax, kmax, P);
        setBoundaryConditionsPressureInDomainOpencl(eargs3D, imax, jmax, kmax, P, Flag);

        copyPressureOpencl(eargsWholeDomain3D, imax, jmax, kmax, P, P_temp);

        sorSolverIterationOpencl(eargs3D, omg, dx, dy, dz, coeff, imax, jmax, kmax, P, P_temp, RS, Flag);
        sorSolverComputeResidualArrayOpencl(
                eargs3D, dx, dy, dz, imax, jmax, kmax, P, RS, Flag,
                openclReductionArrayResidual1, openclReductionArrayNumCells1);

        residual = reduceSumOpenclReal(
                queue, workGroupSize1D, reduceSumOpenclKernelReal, openclReductionArrayResidual1,
                imax*jmax*kmax, openclReductionArrayResidual2, blockSize1D);
        unsigned int numFluidCells = reduceSumOpenclUint(
                queue, workGroupSize1D, reduceSumOpenclKernelUint, openclReductionArrayNumCells1,
                imax*jmax*kmax, openclReductionArrayNumCells2, blockSize1D);
        residual = std::sqrt(residual / Real(numFluidCells));

        it++;
    }

    if (((residual > eps && it == itermax) || std::isnan(residual)) && shallWriteOutput) {
        std::cerr << "\nSOR solver reached maximum number of iterations without converging (res: "
                  << residual << ")." << std::endl;
    }
    if (std::isnan(residual)) {
        std::cerr << "\nResidual in SOR solver is not a number." << std::endl;
        exit(1);
    }
}
