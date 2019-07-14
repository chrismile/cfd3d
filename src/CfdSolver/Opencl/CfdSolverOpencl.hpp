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

#ifndef CFD3D_CFDSOLVEROPENCL_HPP
#define CFD3D_CFDSOLVEROPENCL_HPP

#include <cmath>
#include "CfdSolver/CfdSolver.hpp"
#include "ClInterface.hpp"

class CfdSolverOpencl : public CfdSolver {
public:
    /**
     * @param platformId The ID of the OpenCL platform to use. Which platform corresponds to which ID can be found out
     * with the command line tool 'clinfo'.
     * @param blockSizeX The block size to use for 3D domains in x direction.
     * @param blockSizeY The block size to use for 3D domains in y direction.
     * @param blockSizeZ The block size to use for 3D domains in z direction.
     * @param blockSize1D The block size to use for 1D domains.
     */
    CfdSolverOpencl(int platformId, int blockSizeX, int blockSizeY, int blockSizeZ, int blockSize1D);

    /**
     * Copies the passed initial values of U, V, W, P, T and Flag to the internal representation of the solver.
     * @param scenarioName The name of the scenario as a short string.
     * @param linearSystemSolverType The type of solver to use for solving the Pressure Poisson Equation (PPE).
     * @param Re The Reynolds number used for the simulation.
     * @param Pr The Prandtl number used for the simulation.
     * @param omg The over-relaxation factor of the SOR solver.
     * @param eps The residual value (epsilon) for which the solution of the SOR solver is considered as converged.
     * @param itermax The maximum number of iterations the SOR solver performs until it gives up.
     * @param alpha Donor-cell scheme factor.
     * @param beta Coefficient of thermal expansion.
     * @param dt If tau > 0: The constant time step to use for simulating. Otherwise, dt is overwritten each iteration.
     * @param tau Safety factor \in (0,1] for the maximum time step computation. If tau < 0, the passed value of dt is
     * used as a constant time step.
     * @param GX The gravity in x direction.
     * @param GY The gravity in y direction.
     * @param GZ The gravity in z direction.
     * @param useTemperature Whether the temperature should also be simulated.
     * @param T_h The temperature at boundary cells with the hot temperature flag.
     * @param T_c The temperature at boundary cells with the cold temperature flag.
     * @param imax Number of cells in x direction inside of the domain.
     * @param jmax Number of cells in y direction inside of the domain.
     * @param kmax Number of cells in z direction inside of the domain.
     * @param dx The cell size in x direction.
     * @param dy The cell size in y direction.
     * @param dz The cell size in z direction.
     * @param U The velocities in x direction.
     * @param V The velocities in y direction.
     * @param W The velocities in z direction.
     * @param P The pressure values.
     * @param T The temperature values.
     * @param Flag The flag values (@see Flag.hpp for more information).
     */
    virtual void initialize(
            const std::string &scenarioName, LinearSystemSolverType linearSystemSolverType,
            Real Re, Real Pr, Real omg, Real eps, int itermax, Real alpha, Real beta, Real dt, Real tau,
            Real GX, Real GY, Real GZ, bool useTemperature, Real T_h, Real T_c,
            int imax, int jmax, int kmax, Real dx, Real dy, Real dz,
            Real *U, Real *V, Real *W, Real *P, Real *T, uint32_t *Flag);

    /**
     * The destructor frees the memory of the internal representations of U, V, W, P, T and Flag.
     */
    ~CfdSolverOpencl();


    /**
     * Sets the boundary condition values of U, V, W and T using the Flag array.
     */
    virtual void setBoundaryValues();

    /**
     * Sets special boundary conditions (typically something like inflow) specific to the different scenarios.
     */
    virtual void setBoundaryValuesScenarioSpecific();

    /**
     * Calculates the largest possible (plus some safety margin) time step for the simulation at the current state.
     * @return The time step.
     */
    virtual Real calculateDt();


    /**
     * Updates the temperature values (using an intermediate copy of the temperature from the last iteration).
     */
    virtual void calculateTemperature();

    /**
     * Compute the values in the helper array F, G and H necessary to compute the right-hand-side of the Pressure
     * Poisson equation (PPE).
     */
    virtual void calculateFgh();

    /**
     * Compute the right-hand-side of the Pressure Poisson Equation (PPE).
     */
    virtual void calculateRs();

    /**
     * Execute the SOR solver (successive over-relaxation) for solving the Pressure Poisson Equation (PPE).
     */
    virtual void executeSorSolver();

    /**
     * Updates the values in the arrays U, V and W.
     */
    virtual void calculateUvw();


    /**
     * Copies the values of the internal representations of U, V, W, P and T to the specified arrays.
     * This is necessary when outputting the simulation results at certain time intervals.
     * @param U The velocities in x direction.
     * @param V The velocities in y direction.
     * @param W The velocities in z direction.
     * @param P The pressure values.
     * @param T The temperature values.
     */
    virtual void getDataForOutput(Real *U, Real *V, Real *W, Real *P, Real *T);

private:
    std::string scenarioName;
    LinearSystemSolverType linearSystemSolverType;
    Real Re, Pr, omg, eps, alpha, beta, dt, tau, GX, GY, GZ, T_h, T_c;
    bool useTemperature;
    int itermax;
    int imax, jmax, kmax;
    Real dx, dy, dz;
    cl::Buffer U, V, W, P, P_temp, T, T_temp, F, G, H, RS;
    cl::Buffer Flag;
    int blockSizeX, blockSizeY, blockSizeZ, blockSize1D;
    int platformId;

    cl::Context context;
    std::vector<cl::Device> devices;
    cl::Program computeProgramBoundaryValues, computeProgramSor, computeProgramUvw;
    cl::CommandQueue queue;
    cl::NDRange workGroupSize1D, workGroupSize2D, workGroupSize3D;

    cl::Kernel setLeftRightBoundariesOpenclKernel, setDownUpBoundariesOpenclKernel, setFrontBackBoundariesOpenclKernel,
            setInternalUBoundariesOpenclKernel, setInternalVBoundariesOpenclKernel, setInternalWBoundariesOpenclKernel,
            setInternalTBoundariesOpenclKernel, setDrivenCavityBoundariesOpenclKernel,
            setFlowOverStepBoundariesOpenclKernel, setSingleTowerBoundariesOpenclKernel,
            setMountainBoundariesOpenclKernel;

    cl::Kernel setXYPlanesPressureBoundariesOpenclKernel, setXZPlanesPressureBoundariesOpenclKernel,
            setYZPlanesPressureBoundariesOpenclKernel, setBoundaryConditionsPressureInDomainOpenclKernel,
            copyPressureOpenclKernel, reduceSumOpenclKernelReal, reduceSumOpenclKernelUint,
            sorSolverIterationOpenclKernel, sorSolverComputeResidualArrayOpenclKernel;

    cl::Kernel calculateFghOpenclKernel, setFBoundariesOpenclKernel, setGBoundariesOpenclKernel,
            setHBoundariesOpenclKernel, calculateRsOpenclKernel, calculateMaximumKernel, calculateUvwOpenclKernel,
            calculateTemperatureOpenclKernel;

    // For computing the maximum reduction of the absolute velocities and the sum reduction of the residual.
    cl::Buffer openclReductionArrayU1, openclReductionArrayU2;
    cl::Buffer openclReductionArrayV1, openclReductionArrayV2;
    cl::Buffer openclReductionArrayW1, openclReductionArrayW2;
    cl::Buffer openclReductionArrayResidual1, openclReductionArrayResidual2;
    cl::Buffer openclReductionArrayNumCells1, openclReductionArrayNumCells2;
    cl::Buffer localMemoryReductionReal, localMemoryReductionUint;
};


#endif //CFD3D_CFDSOLVEROPENCL_HPP
