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

#ifndef CFD3D_SORSOLVEROPENCL_HPP
#define CFD3D_SORSOLVEROPENCL_HPP

#include "Defines.hpp"
#include "ClInterface.hpp"

/**
 * Uses an SOR solver to compute the updated pressure values using the pressure poisson equation (PPE).
 */
void sorSolverOpencl(
        Real omg, Real eps, int itermax, LinearSystemSolverType linearSystemSolverType,
        Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        int blockSizeX, int blockSizeY, int blockSizeZ, int blockSize1D,
        cl::CommandQueue &queue, cl::NDRange workGroupSize1D, cl::NDRange workGroupSize2D, cl::NDRange workGroupSize3D,
        cl::Buffer &P, cl::Buffer &P_temp, cl::Buffer &RS, cl::Buffer &Flag,
        cl::Buffer &openclReductionArrayResidual1, cl::Buffer &openclReductionArrayResidual2,
        cl::Buffer &openclReductionArrayNumCells1, cl::Buffer &openclReductionArrayNumCells2,
        cl::Buffer &localMemoryReductionReal, cl::Buffer &localMemoryReductionUint,
        cl::Kernel &setXYPlanesPressureBoundariesOpenclKernel, cl::Kernel &setXZPlanesPressureBoundariesOpenclKernel,
        cl::Kernel &setYZPlanesPressureBoundariesOpenclKernel,
        cl::Kernel &setBoundaryConditionsPressureInDomainOpenclKernel, cl::Kernel &copyPressureOpenclKernel,
        cl::Kernel &reduceSumOpenclKernelReal, cl::Kernel &reduceSumOpenclKernelUint,
        cl::Kernel &sorSolverIterationOpenclKernel, cl::Kernel &sorSolverComputeResidualArrayOpenclKernel);

#endif //CFD3D_SORSOLVEROPENCL_HPP
