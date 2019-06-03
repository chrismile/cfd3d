//
// Created by christoph on 03.06.19.
//

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
void readScenarioConfigurationFromFile(
        const std::string &scenarioFilename, std::string &scenarioName, std::string &geometryName,
        Real &tEnd, Real &dtWrite, Real &xLength, Real &yLength, Real &zLength,
        Real &xOrigin, Real &yOrigin, Real &zOrigin,
        Real &Re, Real &Pr, Real &omg, Real &eps, int &itermax, Real &alpha, Real &beta, Real &dt, Real &tau,
        Real &GX, Real &GY, Real &GZ, bool &useTemperature, Real &T_h, Real &T_c,
        int &imax, int &jmax, int &kmax, Real &dx, Real &dy, Real &dz);

#endif //CFD3D_SCENARIOFILE_HPP
