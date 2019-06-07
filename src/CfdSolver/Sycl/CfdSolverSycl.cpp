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
#include "BoundaryValuesSycl.hpp"
#include "UvwSycl.hpp"
#include "SorSolverSycl.hpp"
#include "CfdSolverSycl.hpp"

CfdSolverSycl::CfdSolverSycl(
        int imax, int jmax, int kmax,
        Real *U, Real *V, Real *W, Real *P, Real *T, uint32_t *Flag)
    : UBuffer(U, cl::sycl::range<1>((imax+1)*(jmax+2)*(kmax+2))),
        VBuffer(V, cl::sycl::range<1>((imax+2)*(jmax+1)*(kmax+2))),
        WBuffer(W, cl::sycl::range<1>((imax+2)*(jmax+2)*(kmax+1))),
        PBuffer(P, cl::sycl::range<1>((imax+2)*(jmax+2)*(kmax+2))),
        TBuffer(T, cl::sycl::range<1>((imax+2)*(jmax+2)*(kmax+2))),
        T_tempBuffer(cl::sycl::range<1>((imax+2)*(jmax+2)*(kmax+2))),
        FBuffer(cl::sycl::range<1>((imax+1)*(jmax+1)*(kmax+1))),
        GBuffer(cl::sycl::range<1>((imax+1)*(jmax+1)*(kmax+1))),
        HBuffer(cl::sycl::range<1>((imax+1)*(jmax+1)*(kmax+1))),
        RSBuffer(cl::sycl::range<1>((imax+1)*(jmax+1)*(kmax+1))),
        FlagBuffer(Flag, cl::sycl::range<1>((imax+2)*(jmax+2)*(kmax+2))) {}

void CfdSolverSycl::initialize(const std::string &scenarioName,
        Real Re, Real Pr, Real omg, Real eps, int itermax, Real alpha, Real beta, Real dt, Real tau,
        Real GX, Real GY, Real GZ, bool useTemperature, Real T_h, Real T_c,
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz,
        Real *U, Real *V, Real *W, Real *P, Real *T, uint32_t *Flag) {
    this->scenarioName = scenarioName;
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

    auto exceptionHandler = [] (cl::sycl::exception_list exceptions) {
        for (std::exception_ptr const& e : exceptions) {
            try {
                std::rethrow_exception(e);
            } catch (cl::sycl::exception const& e) {
                std::cout << "Caught asynchronous SYCL exception:" << std::endl << e.what() << std::endl;
            }
        }
    };

    queue = cl::sycl::queue(cl::sycl::default_selector{}, exceptionHandler);
    cl::sycl::device Device = queue.get_device();
    std::cout << "Device: " << Device.get_info<cl::sycl::info::device::name>() << std::endl;
    std::cout << "Is GPU: " << Device.is_gpu() << std::endl;
}

CfdSolverSycl::~CfdSolverSycl() {
}

void CfdSolverSycl::setBoundaryValues() {
    setBoundaryValuesSycl(queue, T_h, T_c, imax, jmax, kmax,
            UBuffer, VBuffer, WBuffer, TBuffer, FlagBuffer);
}

void CfdSolverSycl::setBoundaryValuesScenarioSpecific() {
    setBoundaryValuesScenarioSpecificSycl(queue, scenarioName, imax, jmax, kmax,
            UBuffer, VBuffer, WBuffer, FlagBuffer);
}

Real CfdSolverSycl::calculateDt() {
    calculateDtSycl(queue, Re, Pr, tau, dt, dx, dy, dz, imax, jmax, kmax,
            UBuffer, VBuffer, WBuffer, useTemperature);
    return dt;
}


void CfdSolverSycl::calculateTemperature() {
    calculateTemperatureSycl(queue, Re, Pr, alpha, dt, dx, dy, dz, imax, jmax, kmax,
            UBuffer, VBuffer, WBuffer, TBuffer, T_tempBuffer, FlagBuffer);
}

void CfdSolverSycl::calculateFgh() {
    calculateFghSycl(queue, Re, GX, GY, GZ, alpha, beta, dt, dx, dy, dz, imax, jmax, kmax,
            UBuffer, VBuffer, WBuffer, TBuffer, FBuffer, GBuffer, HBuffer, FlagBuffer);
}

void CfdSolverSycl::calculateRs() {
    calculateRsSycl(queue, dt, dx, dy, dz, imax, jmax, kmax, FBuffer, GBuffer, HBuffer, RSBuffer);
}


void CfdSolverSycl::executeSorSolver() {
    sorSolverSycl(queue, omg, eps, itermax, dx, dy, dz, imax, jmax, kmax, PBuffer, RSBuffer, FlagBuffer);
}

void CfdSolverSycl::calculateUvw() {
    calculateUvwSycl(queue, dt, dx, dy, dz, imax, jmax, kmax, UBuffer, VBuffer, WBuffer,
            FBuffer, GBuffer, HBuffer, PBuffer, FlagBuffer);
}

void CfdSolverSycl::getDataForOutput(Real *U, Real *V, Real *W, Real *P, Real *T) {
    // Copy the content of U, V, W, P, T in the internal representation to the specified output arrays.
    queue.submit([&](cl::sycl::handler &cgh) {
        ReadAccReal URead = UBuffer.get_access<cl::sycl::access::mode::read_write>(cgh);
        ReadAccReal VRead = VBuffer.get_access<cl::sycl::access::mode::read_write>(cgh);
        ReadAccReal WRead = WBuffer.get_access<cl::sycl::access::mode::read_write>(cgh);
        ReadAccReal PRead = PBuffer.get_access<cl::sycl::access::mode::read_write>(cgh);
        ReadAccReal TRead = TBuffer.get_access<cl::sycl::access::mode::read_write>(cgh);
        cgh.copy(URead, U);
        cgh.copy(VRead, V);
        cgh.copy(WRead, W);
        cgh.copy(PRead, P);
        cgh.copy(TRead, T);
    });

    try {
        queue.wait_and_throw();
    } catch (cl::sycl::exception const& exception) {
        std::cout << "Caught synchronous SYCL exception:n" << std::endl << exception.what() << std::endl;
    }
}
