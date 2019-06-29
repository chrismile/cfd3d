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
#include "SorSolverSycl.hpp"

void sorSolverIterationSycl(
        cl::sycl::queue &queue,
        Real omg, Real dx, Real dy, Real dz, Real coeff, int imax, int jmax, int kmax,
        cl::sycl::buffer<Real, 1> &PBuffer, cl::sycl::buffer<Real, 1> &P_tempBuffer,
        cl::sycl::buffer<Real, 1> &RSBuffer, cl::sycl::buffer<unsigned int, 1> &FlagBuffer,
        Real &residual) {
    // TODO
}

void sorSolverSycl(
        cl::sycl::queue &queue,
        Real omg, Real eps, int itermax,
        Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        cl::sycl::buffer<Real, 1> &PBuffer, cl::sycl::buffer<Real, 1> &P_tempBuffer,
        cl::sycl::buffer<Real, 1> &RSBuffer, cl::sycl::buffer<unsigned int, 1> &FlagBuffer) {
    const Real coeff = omg / (2.0 * (1.0 / (dx*dx) + 1.0 / (dy*dy) + 1.0 / (dz*dz)));
    Real residual = 1e9;
    int it = 0;
    while (it < itermax && residual > eps) {
        queue.submit([&](cl::sycl::handler &cgh) {
            ReadAccReal PRead = PBuffer.get_access<cl::sycl::access::mode::read>(cgh);
            ReadAccReal P_tempWrite = P_tempBuffer.get_access<cl::sycl::access::mode::write>(cgh);
            cgh.copy(PRead, P_tempWrite);
        });
        sorSolverIterationSycl(omg, dx, dy, dz, coeff, imax, jmax, kmax, P, P_temp, RS, Flag, residual);
        it++;
    }

    if (residual > eps && it == itermax) {
        std::cerr << "\nSOR solver reached maximum number of iterations without converging." << std::endl;
    }
}
