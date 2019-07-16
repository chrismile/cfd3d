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

#ifndef CFD3D_UVWOPENCL_HPP
#define CFD3D_UVWOPENCL_HPP

#include "Defines.hpp"
#include "ClInterface.hpp"

/*
 * Determines the maximum time step size. The time step size is restricted according to the CFL theorem.
 */
void calculateDtOpencl(
        Real Re, Real Pr, Real tau,
        Real &dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax, int blockSize1D,
        cl::CommandQueue &queue, cl::NDRange workGroupSize1D, cl::Kernel &calculateMaximumKernel,
        cl::Buffer &U, cl::Buffer &V, cl::Buffer &W,
        cl::Buffer &openclReductionArrayU1, cl::Buffer &openclReductionArrayU2,
        cl::Buffer &openclReductionArrayV1, cl::Buffer &openclReductionArrayV2,
        cl::Buffer &openclReductionArrayW1, cl::Buffer &openclReductionArrayW2,
        bool useTemperature);

#endif //CFD3D_UVWOPENCL_HPP
