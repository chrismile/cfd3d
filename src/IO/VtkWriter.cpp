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

#include <cstdio>
#include <iostream>
#include "VtkWriter.hpp"

bool VtkWriter::initializeWriter(const std::string &filename,
        int imax, int jmax, int kmax, Real dx, Real dy, Real dz, Real xOrigin, Real yOrigin, Real zOrigin) {
    this->filename = filename;
    this->imax = imax;
    this->jmax = jmax;
    this->kmax = kmax;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
    this->xOrigin = xOrigin;
    this->yOrigin = yOrigin;
    this->zOrigin = zOrigin;

    pointData = new Real[(imax+1)*(jmax+1)*(kmax+1)*3];
    cellData = new Real[(imax+1)*(jmax+1)*(kmax+1)];
    cellDataUint = new uint8_t[(imax+1)*(jmax+1)*(kmax+1)];

    return true;
}

VtkWriter::~VtkWriter() {
    if (pointData != nullptr) {
        delete[] pointData;
        delete[] cellData;
        delete[] cellDataUint;
    }
}

void VtkWriter::writeTimestep(
        int timeStepNumber, Real time, Real *U, Real *V, Real *W, Real *P, Real *T, FlagType *Flag) {
    std::string vtkFilename = filename + "." + std::to_string(timeStepNumber) + ".vtk";
    FILE *file = nullptr;
    if (isBinaryVtk) {
        file = fopen(vtkFilename.c_str(), "wb");
    } else {
        file = fopen(vtkFilename.c_str(), "w");
    }
    if (file == nullptr) {
        std::cerr << "Error: Couldn't open file \"" << vtkFilename << "\" for writing." << std::endl;
        exit(1);
    }

    writeVtkHeader(file, imax, jmax, kmax);
    writePointCoordinates(file, imax, jmax, kmax, dx, dy, dz, xOrigin, yOrigin, zOrigin);
    writePointData(file, imax, jmax, kmax, U, V, W, Flag);
    writeCellData(file, imax, jmax, kmax, P, T, Flag);

    fclose(file);
}

void VtkWriter::writeVtkHeader(FILE *file, int imax, int jmax, int kmax) {
    fprintf(file, "# vtk DataFile Version 2.0\n");
    fprintf(file, "Generated by cfd3d\n");
    if (isBinaryVtk) {
        fprintf(file, "BINARY\n");
    } else {
        fprintf(file, "ASCII\n");
    }
    fprintf(file, "DATASET STRUCTURED_GRID\n");
    fprintf(file, "DIMENSIONS %i %i %i\n", imax+1, jmax+1, kmax+1);
}


#define IDXPT(i,j,k,l) (((k)*(jmax+1)*(imax+1) + (j)*(imax+1) + (i))*3 + (l))
#define IDXCELL(i,j,k) ((k)*(jmax+1)*(imax+1) + (j)*(imax+1) + (i))

template <typename T>
void swapEndianness(T *values, int n) {
    uint8_t *byteArray = (uint8_t*)values;
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        uint8_t swappedEntry[sizeof(T)];
        for (int j = 0; j < sizeof(T); j++) {
            swappedEntry[j] = byteArray[i*sizeof(T) + sizeof(T)-j-1];
        }
        values[i] = *((T*)&swappedEntry);
    }
}

void VtkWriter::writePointCoordinates(FILE *file, int imax, int jmax, int kmax,
        Real dx, Real dy, Real dz, Real xOrigin, Real yOrigin, Real zOrigin) {
    fprintf(file, "POINTS %i float\n", (imax+1)*(jmax+1)*(kmax+1));

    if (isBinaryVtk) {
        #pragma omp parallel for
        for (int k = 0; k <= kmax; k++) {
            for (int j = 0; j <= jmax; j++) {
                for (int i = 0; i <= imax; i++) {
                    pointData[IDXPT(i,j,k,0)] = xOrigin + i * dx;
                    pointData[IDXPT(i,j,k,1)] = yOrigin + j * dy;
                    pointData[IDXPT(i,j,k,2)] = zOrigin + k * dz;
                }
            }
        }
        swapEndianness(pointData, (imax+1)*(jmax+1)*(kmax+1)*3);
        fwrite(pointData, sizeof(float), (imax+1)*(jmax+1)*(kmax+1)*3, file);
    } else {
        for (int k = 0; k <= kmax; k++) {
            for (int j = 0; j <= jmax; j++) {
                for (int i = 0; i <= imax; i++) {
                    fprintf(file, "%f %f %f\n", xOrigin + i * dx, yOrigin + j * dy, zOrigin + k * dz);
                }
            }
        }
    }
}

void VtkWriter::writePointData(FILE *file, int imax, int jmax, int kmax, Real *U, Real *V, Real *W, FlagType *Flag) {
    fprintf(file, "POINT_DATA %i\n", (imax+1)*(jmax+1)*(kmax+1));
    fprintf(file, "VECTORS velocity float\n");

    if (isBinaryVtk) {
        #pragma omp parallel for
        for (int k = 0; k <= kmax; k++) {
            for (int j = 0; j <= jmax; j++) {
                for (int i = 0; i <= imax; i++) {
                    Real u = 0, v = 0, w = 0;
                    if (isFluid(Flag[IDXFLAG(i,j,k)])) {
                        u = (U[IDXU(i,j,k)] + U[IDXU(i,j+1,k)] + U[IDXU(i,j,k+1)] + U[IDXU(i,j+1,k+1)]) * Real(0.25);
                        v = (V[IDXV(i,j,k)] + V[IDXV(i+1,j,k)] + V[IDXV(i,j,k+1)] + V[IDXV(i+1,j,k+1)]) * Real(0.25);
                        w = (W[IDXW(i,j,k)] + W[IDXW(i,j+1,k)] + W[IDXW(i+1,j,k)] + W[IDXW(i+1,j+1,k)]) * Real(0.25);
                    }
                    pointData[IDXPT(i,j,k,0)] = u;
                    pointData[IDXPT(i,j,k,1)] = v;
                    pointData[IDXPT(i,j,k,2)] = w;
                }
            }
        }
        swapEndianness(pointData, (imax+1)*(jmax+1)*(kmax+1)*3);
        fwrite(pointData, sizeof(float), (imax+1)*(jmax+1)*(kmax+1)*3, file);
    } else {
        for (int k = 0; k <= kmax; k++) {
            for (int j = 0; j <= jmax; j++) {
                for (int i = 0; i <= imax; i++) {
                    if (isFluid(Flag[IDXFLAG(i,j,k)])) {
                        Real u = (U[IDXU(i,j,k)] + U[IDXU(i,j+1,k)] + U[IDXU(i,j,k+1)] + U[IDXU(i,j+1,k+1)]) * Real(0.25);
                        Real v = (V[IDXV(i,j,k)] + V[IDXV(i+1,j,k)] + V[IDXV(i,j,k+1)] + V[IDXV(i+1,j,k+1)]) * Real(0.25);
                        Real w = (W[IDXW(i,j,k)] + W[IDXW(i,j+1,k)] + W[IDXW(i+1,j,k)] + W[IDXW(i+1,j+1,k)]) * Real(0.25);
                        fprintf(file, "%f %f %f\n", u, v, w);
                    } else {
                        fprintf(file, "0 0 0\n");
                    }
                }
            }
        }
    }
}

void VtkWriter::writeCellData(FILE *file, int imax, int jmax, int kmax, Real *P, Real *T, FlagType *Flag) {
    fprintf(file, "CELL_DATA %i\n", imax*jmax*kmax);
    fprintf(file, "SCALARS pressure float 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");
    if (isBinaryVtk) {
        #pragma omp parallel for
        for (int k = 1; k <= kmax; k++) {
            for (int j = 1; j <= jmax; j++) {
                for (int i = 1; i <= imax; i++) {
                    if (isFluid(Flag[IDXFLAG(i,j,k)])) {
                        cellData[IDXCELL(i,j,k)] = P[IDXP(i,j,k)];
                    } else {
                        cellData[IDXCELL(i,j,k)] = 0;
                    }
                }
            }
        }
        swapEndianness(cellData, imax*jmax*kmax);
        fwrite(cellData, sizeof(float), imax*jmax*kmax, file);
    } else {
        for (int k = 1; k <= kmax; k++) {
            for (int j = 1; j <= jmax; j++) {
                for (int i = 1; i <= imax; i++) {
                    if (isFluid(Flag[IDXFLAG(i,j,k)])) {
                        fprintf(file, "%f\n", P[IDXP(i,j,k)]);
                    } else {
                        fprintf(file, "0\n");
                    }
                }
            }
        }
    }

    fprintf(file, "SCALARS temperature float 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");
    if (isBinaryVtk) {
        #pragma omp parallel for
        for (int k = 1; k <= kmax; k++) {
            for (int j = 1; j <= jmax; j++) {
                for (int i = 1; i <= imax; i++) {
                    if (isFluid(Flag[IDXFLAG(i,j,k)])) {
                        cellData[IDXCELL(i,j,k)] = T[IDXT(i,j,k)];
                    } else {
                        cellData[IDXCELL(i,j,k)] = 0;
                    }
                }
            }
        }
        swapEndianness(cellData, imax*jmax*kmax);
        fwrite(cellData, sizeof(float), imax*jmax*kmax, file);
    } else {
        for (int k = 1; k <= kmax; k++) {
            for (int j = 1; j <= jmax; j++) {
                for (int i = 1; i <= imax; i++) {
                    if (isFluid(Flag[IDXFLAG(i,j,k)])) {
                        fprintf(file, "%f\n", T[IDXT(i,j,k)]);
                    } else {
                        fprintf(file, "0\n");
                    }
                }
            }
        }
    }

    if (isBinaryVtk) {
        fprintf(file, "SCALARS geometry unsigned_char 1\n");
        fprintf(file, "LOOKUP_TABLE default\n");
        for (int k = 1; k <= kmax; k++) {
            for (int j = 1; j <= jmax; j++) {
                for (int i = 1; i <= imax; i++) {
                    cellDataUint[IDXCELL(i,j,k)] = (Flag[IDXFLAG(i,j,k)] & 0x1) == 1;
                }
            }
        }
        fwrite(cellDataUint, sizeof(uint8_t), imax*jmax*kmax, file);
    } else {
        fprintf(file, "SCALARS geometry bit 1\n");
        fprintf(file, "LOOKUP_TABLE default\n");
        for (int k = 1; k <= kmax; k++) {
            for (int j = 1; j <= jmax; j++) {
                for (int i = 1; i <= imax; i++) {
                    fprintf(file, "%u\n", (Flag[IDXFLAG(i, j, k)] & 0x1) == 1);
                }
            }
        }
    }
}
