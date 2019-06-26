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
#include <netcdf.h>
#include "NetCdfWriter.hpp"

bool NetCdfWriter::openFile(const std::string &filename,
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz, Real xOrigin, Real yOrigin, Real zOrigin) {
    this->imax = imax;
    this->jmax = jmax;
    this->kmax = kmax;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
    this->xOrigin = xOrigin;
    this->yOrigin = yOrigin;
    this->zOrigin = zOrigin;

    // Another file still open?
    if (isFileOpen) {
        delete[] centerCellU;
        delete[] centerCellV;
        delete[] centerCellW;
        nc_close(ncid);
    }

    // Open the NetCDF file for reading
    int status = nc_create(filename.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid);
    if (status != 0) {
        std::cerr << "ERROR in loadNetCdfFile: File \"" << filename << "\" couldn't be opened!" << std::endl;
        return false;
    }
    isFileOpen = true;
    writeIndex = 0;

    ncPutAttributeText(NC_GLOBAL, "Conventions", "CF-1.5");
    ncPutAttributeText(NC_GLOBAL, "title", "Fluid simulation in 3D");
    ncPutAttributeText(NC_GLOBAL, "history", "cfd3d");
    ncPutAttributeText(NC_GLOBAL, "institution", "Technische Universitaet Muenchen, Department of Informatics");
    ncPutAttributeText(NC_GLOBAL, "source", "Simulation using incompressible Navier-Stokes equations.");
    ncPutAttributeText(NC_GLOBAL, "references", "https://gitlab.lrz.de/ga42wis/cfd3d");
    ncPutAttributeText(NC_GLOBAL, "comment", "cfd3d was released under the 3-clause BSD license.");


    // Create dimensions
    int timeDim, xDim, yDim, zDim;
    nc_def_dim(ncid, "time", NC_UNLIMITED, &timeDim);
    nc_def_dim(ncid, "x", imax, &xDim);
    nc_def_dim(ncid, "y", jmax, &yDim);
    nc_def_dim(ncid, "z", kmax, &zDim);

    // Define the time and cell center variables.

    nc_def_var(ncid, "x", NC_REAL, 1, &xDim, &xVar);
    nc_def_var(ncid, "y", NC_REAL, 1, &yDim, &yVar);
    nc_def_var(ncid, "z", NC_REAL, 1, &zDim, &zVar);

    nc_def_var(ncid, "time", NC_REAL, 1, &timeDim, &timeVar);
    ncPutAttributeText(timeVar, "long_name", "Time");
    ncPutAttributeText(timeVar, "units", "seconds since simulation start");

    ncPutAttributeText(xVar, "coordinate_type", "Cartesian X");
    ncPutAttributeText(yVar, "coordinate_type", "Cartesian Y");
    ncPutAttributeText(zVar, "coordinate_type", "Cartesian Z");
    ncPutAttributeText(xVar, "axis", "X");
    ncPutAttributeText(yVar, "axis", "Y");
    ncPutAttributeText(zVar, "axis", "Z");
    ncPutAttributeText(timeVar, "axis", "T");

    // Define the domain variables. The fastest changing index is on the right (C/C++ syntax).
    int dimsTimeIndependent3D[] = {xDim, yDim, zDim};
    int dimsTimeDependent3D[] = {timeDim, xDim, yDim, zDim};
    nc_def_var(ncid, "geometry", NC_REAL, 3, dimsTimeIndependent3D, &geometryVar);
    nc_def_var(ncid, "U", NC_REAL, 4, dimsTimeDependent3D, &UVar);
    nc_def_var(ncid, "V", NC_REAL, 4, dimsTimeDependent3D, &VVar);
    nc_def_var(ncid, "W", NC_REAL, 4, dimsTimeDependent3D, &WVar);
    nc_def_var(ncid, "P", NC_REAL, 4, dimsTimeDependent3D, &PVar);
    nc_def_var(ncid, "T", NC_REAL, 4, dimsTimeDependent3D, &TVar);

    // Write the grid cell centers to the x, y and z variables.
    float gridPosition = xOrigin + Real(0.5) * dx;
    for (size_t i = 0; i < (size_t)imax; i++) {
        nc_put_var1_real(ncid, xVar, &i, &gridPosition);
        gridPosition += dx;
    }
    gridPosition = yOrigin + Real(0.5) * dy;
    for (size_t j = 0; j < (size_t)jmax; j++) {
        nc_put_var1_real(ncid, yVar, &j, &gridPosition);
        gridPosition += dy;
    }
    gridPosition = zOrigin + Real(0.5) * dz;
    for (size_t k = 0; k < (size_t)kmax; k++) {
        nc_put_var1_real(ncid, zVar, &k, &gridPosition);
        gridPosition += dz;
    }

    // Intermediate arrays for computing U, V and W at the cell centers.
    centerCellU = new Real[imax*jmax*kmax];
    centerCellV = new Real[imax*jmax*kmax];
    centerCellW = new Real[imax*jmax*kmax];

    return true;
}

void NetCdfWriter::writeTimeDependentVariable3D_Staggered(int ncVar, int jsize, int ksize, Real *values) {
    size_t start[] = {writeIndex, 0, 0, 0};
    size_t count[] = {1, 1, 1, (size_t)kmax};
    for (int i = 1; i <= imax; i++) {
        start[1] = i-1;
        for (int j = 1; j <= jmax; j++) {
            start[2] = j-1;
            nc_put_vara_real(ncid, ncVar, start, count, (values + i*jsize*ksize + j*ksize + 1));
        }
    }
}

void NetCdfWriter::writeTimeDependentVariable3D_Normal(int ncVar, int jsize, int ksize, Real *values) {
    size_t start[] = {writeIndex, 0, 0, 0};
    size_t count[] = {1, 1, 1, (size_t)kmax};
    for (int i = 0; i < imax; i++) {
        start[1] = i;
        for (int j = 0; j < jmax; j++) {
            start[2] = j;
            nc_put_vara_real(ncid, ncVar, start, count, (values + i*jsize*ksize + j*ksize));
        }
    }
}

NetCdfWriter::~NetCdfWriter() {
    if (isFileOpen) {
        delete[] centerCellU;
        delete[] centerCellV;
        delete[] centerCellW;
        nc_close(ncid);
    }
}

void NetCdfWriter::ncPutAttributeText(int varid, const std::string &name, const std::string &value) {
    nc_put_att_text(ncid, varid, name.c_str(), value.size(), value.c_str());
}

void NetCdfWriter::writeTimestep(int timeStepNumber, Real time, Real *U, Real *V, Real *W, Real *P, Real *T,
        FlagType *Flag) {
    if (timeStepNumber == 0) {
        Real *geometryData = new Real[imax*jmax*kmax];
        #pragma omp parallel for
        for (int i = 0; i < imax; i++) {
            for (int j = 0; j < jmax; j++) {
                for (int k = 0; k < kmax; k++) {
                    geometryData[i*jmax*kmax + j*kmax + k] = isFluid(Flag[IDXFLAG(i+1,j+1,k+1)]) ? 1 : 0;
                }
            }
        }
        nc_put_var_real(ncid, geometryVar, geometryData);
        delete[] geometryData;
    }

    #pragma omp parallel for
    for (int i = 0; i < imax; i++) {
        for (int j = 0; j < jmax; j++) {
            for (int k = 0; k < kmax; k++) {
                centerCellU[i*jmax*kmax + j*kmax + k] = (U[IDXU(i+1,j+1,k+1)] + U[IDXU(i,j+1,k+1)]) / Real(2.0);
                centerCellV[i*jmax*kmax + j*kmax + k] = (V[IDXV(i+1,j+1,k+1)] + V[IDXV(i+1,j,k+1)]) / Real(2.0);
                centerCellW[i*jmax*kmax + j*kmax + k] = (W[IDXW(i+1,j+1,k+1)] + W[IDXW(i+1,j+1,k)]) / Real(2.0);
            }
        }
    }

    // Write the new time
    nc_put_var1_real(ncid, timeVar, &writeIndex, &time);
    writeTimeDependentVariable3D_Normal(UVar, jmax, kmax, centerCellU);
    writeTimeDependentVariable3D_Normal(VVar, jmax, kmax, centerCellV);
    writeTimeDependentVariable3D_Normal(WVar, jmax, kmax, centerCellW);
    writeTimeDependentVariable3D_Staggered(PVar, jmax+2, kmax+2, P);
    writeTimeDependentVariable3D_Staggered(TVar, jmax+2, kmax+2, T);
    writeIndex++;

    //nc_sync(ncid);
}
