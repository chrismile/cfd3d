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

#ifndef CFD3D_BOUNDARYVALUESSYCL_HPP
#define CFD3D_BOUNDARYVALUESSYCL_HPP

#include <string>
#include "Defines.hpp"
#include "SyclDefines.hpp"

/**
 * Sets the boundary condition values of U, V, W and T using the Flag array.
 */
void setBoundaryValuesSycl(
        cl::sycl::queue &queue,
        Real T_h, Real T_c,
        int imax, int jmax, int kmax,
        cl::sycl::buffer<Real, 1> &UBuffer, cl::sycl::buffer<Real, 1> &VBuffer,
        cl::sycl::buffer<Real, 1> &WBuffer, cl::sycl::buffer<Real, 1> &TBuffer,
        cl::sycl::buffer<unsigned int, 1> &FlagBuffer);

/**
 * Sets special boundary conditions (typically something like inflow) specific to the different scenarios.
 */
void setBoundaryValuesScenarioSpecificSycl(
        cl::sycl::queue &queue,
        const std::string &scenarioName,
        int imax, int jmax, int kmax,
        cl::sycl::buffer<Real, 1> &UBuffer, cl::sycl::buffer<Real, 1> &VBuffer, cl::sycl::buffer<Real, 1> &WBuffer,
        cl::sycl::buffer<unsigned int, 1> &FlagBuffer);

#endif //CFD3D_BOUNDARYVALUESSYCL_HPP
