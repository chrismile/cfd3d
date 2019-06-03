//
// Created by christoph on 01.06.19.
//

#include <cstring>
#include "BoundaryValuesCpp.hpp"
#include "UvwCpp.hpp"
#include "SorSolverCpp.hpp"
#include "CfdSolverCpp.hpp"

void CfdSolverCpp::initialize(const std::string &scenarioName,
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

    // Create all arrays for the simulation.
    this->U = new Real[(imax+1)*(jmax+2)*(kmax+2)];
    this->V = new Real[(imax+2)*(jmax+1)*(kmax+2)];
    this->W = new Real[(imax+2)*(jmax+2)*(kmax+1)];
    this->P = new Real[(imax+2)*(jmax+2)*(kmax+2)];
    this->T = new Real[(imax+2)*(jmax+2)*(kmax+2)];
    this->T_temp = new Real[(imax+2)*(jmax+2)*(kmax+2)];
    this->F = new Real[(imax+1)*(jmax+1)*(kmax+1)];
    this->G = new Real[(imax+1)*(jmax+1)*(kmax+1)];
    this->H = new Real[(imax+1)*(jmax+1)*(kmax+1)];
    this->RS = new Real[(imax+1)*(jmax+1)*(kmax+1)];
    this->Flag = new FlagType[(imax+2)*(jmax+2)*(kmax+2)];

    // Copy the content of U, V, W, P, T and Flag to the internal representation.
    memcpy(this->U, U, sizeof(Real)*(imax+1)*(jmax+2)*(kmax+2));
    memcpy(this->V, V, sizeof(Real)*(imax+2)*(jmax+1)*(kmax+2));
    memcpy(this->W, W, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+1));
    memcpy(this->P, P, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2));
    memcpy(this->T, T, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2));
    memcpy(this->Flag, Flag, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2));
}

CfdSolverCpp::~CfdSolverCpp() {
    delete[] U;
    delete[] V;
    delete[] W;
    delete[] P;
    delete[] T;
    delete[] T_temp;
    delete[] F;
    delete[] G;
    delete[] H;
    delete[] RS;
    delete[] Flag;
}

void CfdSolverCpp::setBoundaryValues() {
    setBoundaryValuesCpp(T_h, T_c, imax, jmax, kmax, U, V, W, T, Flag);
}

void CfdSolverCpp::setBoundaryValuesScenarioSpecific() {
    setBoundaryValuesScenarioSpecificCpp(scenarioName, imax, jmax, kmax, U, V, W, Flag);
}

Real CfdSolverCpp::calculateDt() {
    calculateDtCpp(Re, Pr, tau, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, useTemperature);
    return dt;
}


void CfdSolverCpp::calculateTemperature() {
    calculateTemperatureCpp(Re, Pr, alpha, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, T, T_temp, Flag);
}

void CfdSolverCpp::calculateFgh() {
    calculateFghCpp(Re, GX, GY, GZ, alpha, beta, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, T, F, G, H, Flag);
}

void CfdSolverCpp::calculateRs() {
    calculateRsCpp(dt, dx, dy, dz, imax, jmax, kmax, F, G, H, RS);
}


void CfdSolverCpp::executeSorSolver() {
    sorSolverCpp(omg, eps, itermax, dx, dy, dz, imax, jmax, kmax, P, RS, Flag);
}

void CfdSolverCpp::calculateUvw() {
    calculateUvwCpp(dt, dx, dy, dz, imax, jmax, kmax, U, V, W, F, G, H, P, Flag);
}

void CfdSolverCpp::getDataForOutput(Real *U, Real *V, Real *W, Real *P, Real *T) {
    // Copy the content of U, V, W, P, T in the internal representation to the specified output arrays.
    memcpy(U, this->U, sizeof(Real)*(imax+1)*(jmax+2)*(kmax+2));
    memcpy(V, this->V, sizeof(Real)*(imax+2)*(jmax+1)*(kmax+2));
    memcpy(W, this->W, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+1));
    memcpy(P, this->P, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2));
    memcpy(T, this->T, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2));
    memcpy(Flag, this->Flag, sizeof(Real)*(imax+2)*(jmax+2)*(kmax+2));
}
