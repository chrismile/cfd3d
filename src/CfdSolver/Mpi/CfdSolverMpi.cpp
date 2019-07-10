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
#include "BoundaryValuesMpi.hpp"
#include "UvwMpi.hpp"
#include "SorSolverMpi.hpp"
#include "CfdSolverMpi.hpp"
#include "MpiHelpers.hpp"
#include "DefinesMpi.hpp"

CfdSolverMpi::CfdSolverMpi(
        int il, int iu, int jl, int ju, int kl, int ku,
        int myrank, int rankL, int rankR, int rankD, int rankU, int rankB, int rankF) {
    this->il = il;
    this->iu = iu;
    this->jl = jl;
    this->ju = ju;
    this->kl = kl;
    this->ku = ku;
    this->myrank = myrank;
    this->rankL = rankL;
    this->rankR = rankR;
    this->rankD = rankD;
    this->rankU = rankU;
    this->rankB = rankB;
    this->rankF = rankF;
}

void CfdSolverMpi::initialize(const std::string &scenarioName,
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
    this->U = new Real[(iu - il + 4)*(ju - jl + 3)*(ku - kl + 3)];
    this->V = new Real[(iu - il + 3)*(ju - jl + 4)*(ku - kl + 3)];
    this->W = new Real[(iu - il + 3)*(ju - jl + 3)*(ku - kl + 4)];
    this->P = new Real[(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3)];
    this->P_temp = new Real[(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3)];
    this->T = new Real[(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3)];
    this->T_temp = new Real[(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3)];
    this->F = new Real[(iu - il + 4)*(ju - jl + 3)*(ku - kl + 3)];
    this->G = new Real[(iu - il + 3)*(ju - jl + 4)*(ku - kl + 3)];
    this->H = new Real[(iu - il + 3)*(ju - jl + 3)*(ku - kl + 4)];
    this->RS = new Real[(iu - il + 1)*(ju - jl + 1)*(ku - kl + 1)];
    this->Flag = new FlagType[(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3)];

    memset(this->U, 0, sizeof(Real)*(iu - il + 4)*(ju - jl + 3)*(ku - kl + 3));
    memset(this->V, 0, sizeof(Real)*(iu - il + 3)*(ju - jl + 4)*(ku - kl + 3));
    memset(this->W, 0, sizeof(Real)*(iu - il + 3)*(ju - jl + 3)*(ku - kl + 4));
    memset(this->P, 0, sizeof(Real)*(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3));
    memset(this->P_temp, 0, sizeof(Real)*(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3));
    memset(this->T, 0, sizeof(Real)*(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3));
    memset(this->T_temp, 0, sizeof(Real)*(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3));
    memset(this->F, 0, sizeof(Real)*(iu - il + 4)*(ju - jl + 3)*(ku - kl + 3));
    memset(this->G, 0, sizeof(Real)*(iu - il + 3)*(ju - jl + 4)*(ku - kl + 3));
    memset(this->H, 0, sizeof(Real)*(iu - il + 3)*(ju - jl + 3)*(ku - kl + 4));
    memset(this->RS, 0, sizeof(Real)*(iu - il + 1)*(ju - jl + 1)*(ku - kl + 1));
    memset(this->Flag, 0, sizeof(FlagType)*(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3));

    // Copy the content of U, V, W, P, T and Flag to the internal representation.
    memcpy(this->U, U, sizeof(Real)*(iu - il + 4)*(ju - jl + 3)*(ku - kl + 3));
    memcpy(this->V, V, sizeof(Real)*(iu - il + 3)*(ju - jl + 4)*(ku - kl + 3));
    memcpy(this->W, W, sizeof(Real)*(iu - il + 3)*(ju - jl + 3)*(ku - kl + 4));
    memcpy(this->P, P, sizeof(Real)*(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3));
    memcpy(this->T, T, sizeof(Real)*(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3));
    memcpy(this->Flag, Flag, sizeof(unsigned int)*(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3));

    int maxMpiBufferSize = std::max(
            std::max((iu - il + 2) * (ju - jl + 2),
                    (iu - il + 2) * (ku - kl + 2)),
            (ju - jl + 2) * (ku - kl + 2));
    bufSend = new Real[maxMpiBufferSize];
    bufRecv = new Real[maxMpiBufferSize];
}

CfdSolverMpi::~CfdSolverMpi() {
    delete[] U;
    delete[] V;
    delete[] W;
    delete[] P;
    delete[] P_temp;
    delete[] T;
    delete[] T_temp;
    delete[] F;
    delete[] G;
    delete[] H;
    delete[] RS;
    delete[] Flag;

    delete[] bufSend;
    delete[] bufRecv;
}

void CfdSolverMpi::setBoundaryValues() {
    setBoundaryValuesMpi(T_h, T_c, imax, jmax, kmax, il, iu, jl, ju, kl, ku, U, V, W, T, Flag);
}

void CfdSolverMpi::setBoundaryValuesScenarioSpecific() {
    setBoundaryValuesScenarioSpecificMpi(scenarioName, imax, jmax, kmax, il, iu, jl, ju, kl, ku, U, V, W, Flag);
}

Real CfdSolverMpi::calculateDt() {
    calculateDtMpi(Re, Pr, tau, dt, dx, dy, dz, imax, jmax, kmax, il, iu, jl, ju, kl, ku, U, V, W, useTemperature);
    return dt;
}


void CfdSolverMpi::calculateTemperature() {
    Real *temp = T;
    T = T_temp;
    T_temp = temp;
    calculateTemperatureMpi(
            Re, Pr, alpha, dt, dx, dy, dz, imax, jmax, kmax, il, iu, jl, ju, kl, ku,
            rankL, rankR, rankD, rankU, rankB, rankF, bufSend, bufRecv, U, V, W, T, T_temp, Flag);
}

void CfdSolverMpi::calculateFgh() {
    calculateFghMpi(Re, GX, GY, GZ, alpha, beta, dt, dx, dy, dz, imax, jmax, kmax, il, iu, jl, ju, kl, ku,
            U, V, W, T, F, G, H, Flag);
}

void CfdSolverMpi::calculateRs() {
    calculateRsMpi(dt, dx, dy, dz, imax, jmax, kmax, il, iu, jl, ju, kl, ku, F, G, H, RS);
}


void CfdSolverMpi::executeSorSolver() {
    sorSolverMpi(
            myrank, omg, eps, itermax, linearSystemSolverType, dx, dy, dz, imax, jmax, kmax, il, iu, jl, ju, kl, ku,
            rankL, rankR, rankD, rankU, rankB, rankF, bufSend, bufRecv, P, P_temp, RS, Flag);
}

void CfdSolverMpi::calculateUvw() {
    calculateUvwMpi(
            dt, dx, dy, dz, imax, jmax, kmax, il, iu, jl, ju, kl, ku,
            rankL, rankR, rankD, rankU, rankB, rankF, bufSend, bufRecv, U, V, W, F, G, H, P, Flag);
}

void CfdSolverMpi::getDataForOutput(Real *U, Real *V, Real *W, Real *P, Real *T) {
    // Copy the content of U, V, W, P, T in the internal representation to the specified output arrays.
    memcpy(U, this->U, sizeof(Real)*(iu - il + 4)*(ju - jl + 3)*(ku - kl + 3));
    memcpy(V, this->V, sizeof(Real)*(iu - il + 3)*(ju - jl + 4)*(ku - kl + 3));
    memcpy(W, this->W, sizeof(Real)*(iu - il + 3)*(ju - jl + 3)*(ku - kl + 4));
    memcpy(P, this->P, sizeof(Real)*(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3));
    memcpy(T, this->T, sizeof(Real)*(iu - il + 3)*(ju - jl + 3)*(ku - kl + 3));
}
