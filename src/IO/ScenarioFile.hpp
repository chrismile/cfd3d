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

#ifndef CFD3D_SCENARIOFILE_HPP
#define CFD3D_SCENARIOFILE_HPP

#include <string>
#include "Defines.hpp"

/**
 * Loads a scenario configuration file in the .dat format.
 * @param scenarioFilename The file name of the scenario configuration file.
 * @param scenarioName The name of the scenario as a short string.
 * @param geometryName The name of the geometry file to use (or none if no geometry file should be loaded).
 * @param tEnd The end time of the simulation (i.e. at which time to stop).
 * @param dtWrite The time intervals in which the intermediate simulation results should be written to a file.
 * @param xLength The length of the interior of the simulation domain in x direction.
 * @param yLength The length of the interior of the simulation domain in y direction.
 * @param zLength The length of the interior of the simulation domain in z direction.
 * @param xOrigin The origin of the interior domain in world space (x coordinate).
 * @param yOrigin The origin of the interior domain in world space (y coordinate).
 * @param zOrigin The origin of the interior domain in world space (z coordinate).
 * @param UI The initialization value for U.
 * @param VI The initialization value for V.
 * @param WI The initialization value for W.
 * @param PI The initialization value for P.
 * @param TI The initialization value for T.
 * @param GX The gravity in x direction.
 * @param GY The gravity in y direction.
 * @param GZ The gravity in z direction.
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
void readScenarioConfigurationFromFile(
        const std::string &scenarioFilename, std::string &scenarioName, std::string &geometryName,
        Real &tEnd, Real &dtWrite, Real &xLength, Real &yLength, Real &zLength,
        Real &xOrigin, Real &yOrigin, Real &zOrigin,
        Real &UI, Real &VI, Real &WI, Real &PI, Real &TI, Real &GX, Real &GY, Real &GZ,
        Real &Re, Real &Pr, Real &omg, Real &eps, int &itermax, Real &alpha, Real &beta, Real &dt, Real &tau,
        bool &useTemperature, Real &T_h, Real &T_c,
        int &imax, int &jmax, int &kmax, Real &dx, Real &dy, Real &dz);

#endif //CFD3D_SCENARIOFILE_HPP
