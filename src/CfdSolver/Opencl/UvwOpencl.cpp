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

#include <algorithm>
#include <iostream>
#include "UvwOpencl.hpp"

void calculateDtOpencl(
        Real Re, Real Pr, Real tau,
        Real &dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax, int blockSize1D,
        cl::CommandQueue &queue, cl::NDRange workGroupSize1D, cl::Kernel &calculateMaximumKernel,
        cl::Buffer &U, cl::Buffer &V, cl::Buffer &W,
        cl::Buffer &openclReductionArrayU1, cl::Buffer &openclReductionArrayU2,
        cl::Buffer &openclReductionArrayV1, cl::Buffer &openclReductionArrayV2,
        cl::Buffer &openclReductionArrayW1, cl::Buffer &openclReductionArrayW2,
        cl::Buffer &localMemoryReductionReal,
        bool useTemperature) {
    Real uMaxAbs = Real(0.0), vMaxAbs = Real(0.0), wMaxAbs = Real(0.0);

    cl::make_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int> calculateMaximum(calculateMaximumKernel);

    cl::Buffer U_reductionInput = U;
    cl::Buffer V_reductionInput = V;
    cl::Buffer W_reductionInput = W;
    cl::Buffer U_reductionOutput = openclReductionArrayU1;
    cl::Buffer V_reductionOutput = openclReductionArrayV1;
    cl::Buffer W_reductionOutput = openclReductionArrayW1;

    int numberOfBlocksI = (imax+1)*(jmax+2)*(kmax+2);
    int numberOfBlocksJ = (imax+2)*(jmax+1)*(kmax+2);
    int numberOfBlocksK = (imax+2)*(jmax+2)*(kmax+1);
    int inputSizeI;
    int inputSizeJ;
    int inputSizeK;
    bool finished = false;

    int iteration = 0;
    while (!finished) {
        inputSizeI = numberOfBlocksI;
        inputSizeJ = numberOfBlocksJ;
        inputSizeK = numberOfBlocksK;
        numberOfBlocksI = iceil(numberOfBlocksI, blockSize1D*2);
        numberOfBlocksJ = iceil(numberOfBlocksJ, blockSize1D*2);
        numberOfBlocksK = iceil(numberOfBlocksK, blockSize1D*2);

        if (inputSizeI != 1) {
            cl::EnqueueArgs eargs(queue, cl::NullRange, cl::NDRange(numberOfBlocksI*blockSize1D*2), workGroupSize1D);
            calculateMaximum(
                    eargs, U_reductionInput, U_reductionOutput, localMemoryReductionReal, inputSizeI);
            if (iteration % 2 == 0) {
                U_reductionInput = openclReductionArrayU1;
                U_reductionOutput = openclReductionArrayU2;
            } else {
                U_reductionInput = openclReductionArrayU2;
                U_reductionOutput = openclReductionArrayU1;
            }
        }
        if (inputSizeJ != 1) {
            cl::EnqueueArgs eargs(queue, cl::NullRange, cl::NDRange(numberOfBlocksJ*blockSize1D*2), workGroupSize1D);
            calculateMaximum(
                    eargs, V_reductionInput, V_reductionOutput, localMemoryReductionReal, inputSizeJ);
            if (iteration % 2 == 0) {
                V_reductionInput = openclReductionArrayV1;
                V_reductionOutput = openclReductionArrayV2;
            } else {
                V_reductionInput = openclReductionArrayV2;
                V_reductionOutput = openclReductionArrayV1;
            }
        }
        if (inputSizeK != 1) {
            cl::EnqueueArgs eargs(queue, cl::NullRange, cl::NDRange(numberOfBlocksK*blockSize1D*2), workGroupSize1D);
            calculateMaximum(
                    eargs, W_reductionInput, W_reductionOutput, localMemoryReductionReal, inputSizeK);
            if (iteration % 2 == 0) {
                W_reductionInput = openclReductionArrayW1;
                W_reductionOutput = openclReductionArrayW2;
            } else {
                W_reductionInput = openclReductionArrayW2;
                W_reductionOutput = openclReductionArrayW1;
            }
        }


        if (numberOfBlocksI == 1 && numberOfBlocksJ == 1 && numberOfBlocksK == 1) {
            finished = true;
        }
        iteration++;
    }

    queue.enqueueReadBuffer(U_reductionInput, CL_FALSE, 0, sizeof(Real), (void*)&uMaxAbs);
    queue.enqueueReadBuffer(V_reductionInput, CL_FALSE, 0, sizeof(Real), (void*)&vMaxAbs);
    queue.enqueueReadBuffer(W_reductionInput, CL_FALSE, 0, sizeof(Real), (void*)&wMaxAbs);
    queue.finish();

    if (tau < Real(0.0)) {
        // Constant time step manually specified in configuration file. Check for stability.
        assert(2 / Re * dt < dx * dx * dy * dy * dz * dz / (dx * dx + dy * dy + dz * dz));
        assert(uMaxAbs * dt < dx);
        assert(vMaxAbs * dt < dy);
        assert(wMaxAbs * dt < dz);
        if (useTemperature){
            assert(dt < (Re*Pr/2)*(1/((1/(dx*dx))+1/(dy*dy)+1/(dz*dz))));
        }
        return;
    }

    dt = std::min(dx / uMaxAbs, dy / vMaxAbs);
    dt = std::min(dt, dz / wMaxAbs);
    dt = std::min(dt, (Re / Real(2.0)) * (Real(1.0) / (Real(1.0) / (dx*dx) + Real(1.0) / (dy*dy)
            + Real(1.0) / (dz*dz))));
    if (useTemperature){
        dt = std::min(dt, (Re * Pr / Real(2.0)) * (Real(1.0) / (Real(1.0) / (dx*dx) + Real(1.0) / (dy*dy)
                + Real(1.0) / (dz*dz))));
    }
    dt = tau * dt;
}

