/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2019, Christoph Neuhauser
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

#ifndef CFD3D_VTKWRITER_HPP
#define CFD3D_VTKWRITER_HPP

#include "OutputFileWriter.hpp"

class VtkWriter : public OutputFileWriter {
public:
    VtkWriter(int nproc, int myrank, bool isBinaryVtk = true) : nproc(nproc), myrank(myrank), isBinaryVtk(isBinaryVtk)
    {}
    virtual ~VtkWriter();

    /**
     * @return The file ending of the format.
     */
    virtual std::string getOutputFormatEnding() { return ".vtk"; }

    /**
     * Sets data for use with the MPI solver.
     */
    virtual void setMpiData(int il, int iu, int jl, int ju, int kl, int ku);

    /**
     * Intializes the file writer.
     * @param filename The file name of the file to write to.
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
    virtual bool initializeWriter(const std::string &filename,
            int imax, int jmax, int kmax, Real dx, Real dy, Real dz, Real xOrigin, Real yOrigin, Real zOrigin);
    /**
     * Writes the data of the current time step to a file.
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
    virtual void writeTimestep(
            int timeStepNumber, Real time, Real *U, Real *V, Real *W, Real *P, Real *T, FlagType *Flag);

private:
    void writeVtkHeader(FILE *file);
    void writePointCoordinates(FILE *file,
            Real dx, Real dy, Real dz, Real xOrigin, Real yOrigin, Real zOrigin);
    void writePointData(FILE *file, Real *U, Real *V, Real *W, FlagType *Flag);
    void writeCellData(FILE *file, Real *P, Real *T, FlagType *Flag);

    bool isBinaryVtk;
    std::string filename;

    bool isMpiMode = false;
    int nproc, myrank;
    int imax, jmax, kmax;
    int il, iu, jl, ju, kl, ku;
    Real dx, dy, dz, xOrigin, yOrigin, zOrigin;

    // For binary mode. (Note: ParaView only supports float, and not double.)
    float *pointData = nullptr;
    float *cellData = nullptr;
    uint8_t *cellDataUint = nullptr;
};


#endif //CFD3D_VTKWRITER_HPP
