//
// Created by christoph on 01.06.19.
//

#ifndef CFD3D_NETCDFWRITER_HPP
#define CFD3D_NETCDFWRITER_HPP

#include <string>
#include "CfdSolver/Flag.hpp"

class NetCdfWriter {
public:
    ~NetCdfWriter();

    /**
     * Opens a NetCDF file for writing. TODO: More about format.
     * @param filename The file name of the NetCDF file to write to.
     * @param imax Number of cells in x direction inside of the domain.
     * @param jmax Number of cells in y direction inside of the domain.
     * @param kmax Number of cells in z direction inside of the domain.
     * @param dx The cell size in x direction.
     * @param dy The cell size in y direction.
     * @param dz The cell size in z direction.
     * @param xOrigin The origin of the interior domain in world space (x coordinate).
     * @param yOrigin The origin of the interior domain in world space (y coordinate).
     * @param zOrigin The origin of the interior domain in world space (z coordinate).
     * @return true if the file could be opened for writing successfully.
     */
    bool openFile(const std::string &filename,
            int imax, int jmax, int kmax, Real dx, Real dy, Real dz, Real xOrigin, Real yOrigin, Real zOrigin);
    /**
     * TODO: More about format.
     * @param timeStepNumber The time step number (i.e., 0 for the initial state, 1 after the first iteration of the
     * simulation, etc.).
     * @param time The current time of the simulation.
     * @param U The velocities in x direction.
     * @param V The velocities in y direction.
     * @param W The velocities in z direction.
     * @param P The pressure values.
     * @param T The temperature values.
     * @param Flag The flag values (@see Flag.hpp for more information).
     */
    void writeTimestep(int timeStepNumber, Real time, Real *U, Real *V, Real *W, Real *P, Real *T, FlagType *Flag);

private:
    void writeTimeDependentVariable3D_Staggered(int timeStepNumber, int ncVar, int jsize, int ksize, Real *values);
    void writeTimeDependentVariable3D_Normal(int timeStepNumber, int ncVar, int jsize, int ksize, Real *values);
    void ncPutAttributeText(int varid, const std::string &name, const std::string &value);
    int imax, jmax, kmax;
    Real dx, dy, dz, xOrigin, yOrigin, zOrigin;
    Real *centerCellU;
    Real *centerCellV;
    Real *centerCellW;

    int ncid;

    // NetCDF variables
    int timeVar, xVar, yVar, zVar, geometryVar, UVar, VVar, WVar, PVar, TVar;
};


#endif //CFD3D_NETCDFWRITER_HPP
