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

#include "UvwSycl.hpp"

void calculateFghSycl(
        cl::sycl::queue &queue,
        Real Re, Real GX, Real GY, Real GZ, Real alpha, Real beta,
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        cl::sycl::buffer<Real, 1> &U, cl::sycl::buffer<Real, 1> &V,
        cl::sycl::buffer<Real, 1> &W, cl::sycl::buffer<Real, 1> &T,
        cl::sycl::buffer<Real, 1> &F, cl::sycl::buffer<Real, 1> &G,
        cl::sycl::buffer<Real, 1> &H, cl::sycl::buffer<unsigned int, 1> &Flag) {
    // TODO
}

void calculateRsSycl(
        cl::sycl::queue &queue,
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        cl::sycl::buffer<Real, 1> &F, cl::sycl::buffer<Real, 1> &G,
        cl::sycl::buffer<Real, 1> &H, cl::sycl::buffer<Real, 1> &RS) {
    // TODO
}

void calculateDtSycl(
        cl::sycl::queue &queue,
        Real Re, Real Pr, Real tau,
        Real &dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        cl::sycl::buffer<Real, 1> &U, cl::sycl::buffer<Real, 1> &V, cl::sycl::buffer<Real, 1> &W,
        bool useTemperature) {
    // TODO
}

void calculateUvwSycl(
        cl::sycl::queue &queue,
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        cl::sycl::buffer<Real, 1> &U, cl::sycl::buffer<Real, 1> &V,
        cl::sycl::buffer<Real, 1> &W, cl::sycl::buffer<Real, 1> &F,
        cl::sycl::buffer<Real, 1> &G, cl::sycl::buffer<Real, 1> &H,
        cl::sycl::buffer<Real, 1> &P, cl::sycl::buffer<unsigned int, 1> &Flag) {
    // TODO
}

void calculateTemperatureSycl(
        cl::sycl::queue &queue,
        Real Re, Real Pr, Real alpha,
        Real dt, Real dx, Real dy, Real dz,
        int imax, int jmax, int kmax,
        cl::sycl::buffer<Real, 1> &U, cl::sycl::buffer<Real, 1> &V,
        cl::sycl::buffer<Real, 1> &W, cl::sycl::buffer<Real, 1> &T,
        cl::sycl::buffer<Real, 1> &T_temp, cl::sycl::buffer<unsigned int, 1> &Flag) {
    // TODO
}
