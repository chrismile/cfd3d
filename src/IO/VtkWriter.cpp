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
#include "CfdSolver/Mpi/DefinesMpi.hpp"

void VtkWriter::setMpiData(int il, int iu, int jl, int ju, int kl, int ku) {
    this->il = il;
    this->iu = iu;
    this->jl = jl;
    this->ju = ju;
    this->kl = kl;
    this->ku = ku;
    isMpiMode = true;
}

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

    if (!isMpiMode) {
        il = jl = kl = 1;
        iu = imax;
        ju = jmax;
        ku = kmax;
    }
    pointData = new float[(iu-il+2)*(ju-jl+2)*(ku-kl+2)*3];
    cellData = new float[(iu-il+2)*(ju-jl+2)*(ku-kl+2)];
    cellDataUint = new uint8_t[(iu-il+2)*(ju-jl+2)*(ku-kl+2)];

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
    if (nproc != 1) {
        // Each process outputs its own file.
        vtkFilename = filename + "." + std::to_string(myrank) + "." + std::to_string(timeStepNumber) + ".vtk";
    }
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

    writeVtkHeader(file);
    writePointCoordinates(file, dx, dy, dz, xOrigin, yOrigin, zOrigin);
    writePointData(file, U, V, W, Flag);
    writeCellData(file, P, T, Flag);

    fclose(file);
}

void VtkWriter::writeVtkHeader(FILE *file) {
    fprintf(file, "# vtk DataFile Version 2.0\n");
    fprintf(file, "Generated by cfd3d\n");
    if (isBinaryVtk) {
        fprintf(file, "BINARY\n");
    } else {
        fprintf(file, "ASCII\n");
    }
    fprintf(file, "DATASET STRUCTURED_GRID\n");
    fprintf(file, "DIMENSIONS %i %i %i\n", iu-il+2, ju-jl+2, ku-kl+2);
}


//#define IDXPT(i,j,k,l) (((k)*(jmax+1)*(imax+1) + (j)*(imax+1) + (i))*3 + (l))
//#define IDXCELL(i,j,k) (((k)-1)*jmax*imax + ((j)-1)*imax + ((i)-1))
#define IDXPT(i,j,k,l) ((((k) - (kl-1))*(ju-jl+2)*(iu-il+2) + ((j) - (jl-1))*(iu-il+2) + ((i) - (il-1)))*3 + (l))
#define IDXCELL(i,j,k) (((k) - kl)*(ju-jl+1)*(iu-il+1) + ((j) - jl)*(iu-il+1) + ((i) - il))

template <typename T>
void swapEndianness(T *values, int n) {
    auto *byteArray = (uint8_t*)values;
    #pragma omp parallel for shared(byteArray, values, n) default(none)
    for (int i = 0; i < n; i++) {
        uint8_t swappedEntry[sizeof(T)];
        for (size_t j = 0; j < sizeof(T); j++) {
            swappedEntry[j] = byteArray[i * sizeof(T) + sizeof(T) - j - 1];
        }
        values[i] = *((T*)&swappedEntry);
    }
}

void VtkWriter::writePointCoordinates(FILE *file,
        Real dx, Real dy, Real dz, Real xOrigin, Real yOrigin, Real zOrigin) {
    fprintf(file, "POINTS %i float\n", (iu-il+2)*(ju-jl+2)*(ku-kl+2));

    if (isBinaryVtk) {
        #pragma omp parallel for shared(dx, dy, dz, xOrigin, yOrigin, zOrigin) default(none)
        for (int k = kl-1; k <= ku; k++) {
            for (int j = jl-1; j <= ju; j++) {
                for (int i = il-1; i <= iu; i++) {
                    pointData[IDXPT(i,j,k,0)] = xOrigin + i * dx;
                    pointData[IDXPT(i,j,k,1)] = yOrigin + j * dy;
                    pointData[IDXPT(i,j,k,2)] = zOrigin + k * dz;
                }
            }
        }
        swapEndianness(pointData, (iu-il+2)*(ju-jl+2)*(ku-kl+2)*3);
        fwrite(pointData, sizeof(float), (iu-il+2)*(ju-jl+2)*(ku-kl+2)*3, file);
    } else {
        for (int k = kl-1; k <= ku; k++) {
            for (int j = jl-1; j <= ju; j++) {
                for (int i = il-1; i <= iu; i++) {
                    fprintf(file, "%f %f %f\n", xOrigin + i * dx, yOrigin + j * dy, zOrigin + k * dz);
                }
            }
        }
    }
}

void VtkWriter::writePointData(FILE *file, Real *U, Real *V, Real *W, FlagType *Flag) {
    fprintf(file, "POINT_DATA %i\n", (iu-il+2)*(ju-jl+2)*(ku-kl+2));
    fprintf(file, "VECTORS velocity float\n");

    if (isBinaryVtk) {
        if (isMpiMode) {
            for (int k = kl-1; k <= ku; k++) {
                for (int j = jl-1; j <= ju; j++) {
                    for (int i = il-1; i <= iu; i++) {
                        Real u = 0, v = 0, w = 0;
                        if (isFluid(Flag[IDXFLAG_MPI(i,j,k)])) {
                            u = (U[IDXU_MPI(i,j,k)] + U[IDXU_MPI(i,j+1,k)] + U[IDXU_MPI(i,j,k+1)]
                                    + U[IDXU_MPI(i,j+1,k+1)]) * Real(0.25);
                            v = (V[IDXV_MPI(i,j,k)] + V[IDXV_MPI(i+1,j,k)] + V[IDXV_MPI(i,j,k+1)]
                                    + V[IDXV_MPI(i+1,j,k+1)]) * Real(0.25);
                            w = (W[IDXW_MPI(i,j,k)] + W[IDXW_MPI(i,j+1,k)] + W[IDXW_MPI(i+1,j,k)]
                                    + W[IDXW_MPI(i+1,j+1,k)]) * Real(0.25);
                        }
                        pointData[IDXPT(i,j,k,0)] = u;
                        pointData[IDXPT(i,j,k,1)] = v;
                        pointData[IDXPT(i,j,k,2)] = w;
                    }
                }
            }
        } else {
            #pragma omp parallel for shared(U, V, W, Flag) default(none)
            for (int k = kl-1; k <= ku; k++) {
                for (int j = jl-1; j <= ju; j++) {
                    for (int i = il-1; i <= iu; i++) {
                        Real u = 0, v = 0, w = 0;
                        if (isFluid(Flag[IDXFLAG_NORMAL(i,j,k)])) {
                            u = (U[IDXU_NORMAL(i,j,k)] + U[IDXU_NORMAL(i,j+1,k)] + U[IDXU_NORMAL(i,j,k+1)]
                                    + U[IDXU_NORMAL(i,j+1,k+1)]) * Real(0.25);
                            v = (V[IDXV_NORMAL(i,j,k)] + V[IDXV_NORMAL(i+1,j,k)] + V[IDXV_NORMAL(i,j,k+1)]
                                    + V[IDXV_NORMAL(i+1,j,k+1)]) * Real(0.25);
                            w = (W[IDXW_NORMAL(i,j,k)] + W[IDXW_NORMAL(i,j+1,k)] + W[IDXW_NORMAL(i+1,j,k)]
                                    + W[IDXW_NORMAL(i+1,j+1,k)]) * Real(0.25);
                        }
                        pointData[IDXPT(i,j,k,0)] = u;
                        pointData[IDXPT(i,j,k,1)] = v;
                        pointData[IDXPT(i,j,k,2)] = w;
                    }
                }
            }
        }
        swapEndianness(pointData, (iu-il+2)*(ju-jl+2)*(ku-kl+2)*3);
        fwrite(pointData, sizeof(float), (iu-il+2)*(ju-jl+2)*(ku-kl+2)*3, file);
    } else {
        if (isMpiMode) {
            for (int k = kl-1; k <= ku; k++) {
                for (int j = jl-1; j <= ju; j++) {
                    for (int i = il-1; i <= iu; i++) {
                        if (isFluid(Flag[IDXFLAG_MPI(i,j,k)])) {
                            Real u = (U[IDXU_MPI(i,j,k)] + U[IDXU_MPI(i,j+1,k)] + U[IDXU_MPI(i,j,k+1)]
                                    + U[IDXU_MPI(i,j+1,k+1)]) * Real(0.25);
                            Real v = (V[IDXV_MPI(i,j,k)] + V[IDXV_MPI(i+1,j,k)] + V[IDXV_MPI(i,j,k+1)]
                                    + V[IDXV_MPI(i+1,j,k+1)]) * Real(0.25);
                            Real w = (W[IDXW_MPI(i,j,k)] + W[IDXW_MPI(i,j+1,k)] + W[IDXW_MPI(i+1,j,k)]
                                    + W[IDXW_MPI(i+1,j+1,k)]) * Real(0.25);
                            fprintf(file, "%f %f %f\n", u, v, w);
                        } else {
                            fprintf(file, "0 0 0\n");
                        }
                    }
                }
            }
        } else {
            for (int k = kl-1; k <= ku; k++) {
                for (int j = jl-1; j <= ju; j++) {
                    for (int i = il-1; i <= iu; i++) {
                        if (isFluid(Flag[IDXFLAG_NORMAL(i,j,k)])) {
                            Real u = (U[IDXU_NORMAL(i,j,k)] + U[IDXU_NORMAL(i,j+1,k)] + U[IDXU_NORMAL(i,j,k+1)]
                                    + U[IDXU_NORMAL(i,j+1,k+1)]) * Real(0.25);
                            Real v = (V[IDXV_NORMAL(i,j,k)] + V[IDXV_NORMAL(i+1,j,k)] + V[IDXV_NORMAL(i,j,k+1)]
                                    + V[IDXV_NORMAL(i+1,j,k+1)]) * Real(0.25);
                            Real w = (W[IDXW_NORMAL(i,j,k)] + W[IDXW_NORMAL(i,j+1,k)] + W[IDXW_NORMAL(i+1,j,k)]
                                    + W[IDXW_NORMAL(i+1,j+1,k)]) * Real(0.25);
                            fprintf(file, "%f %f %f\n", u, v, w);
                        } else {
                            fprintf(file, "0 0 0\n");
                        }
                    }
                }
            }
        }
    }
}

void VtkWriter::writeCellData(FILE *file, Real *P, Real *T, FlagType *Flag) {
    if (isMpiMode) {
        fprintf(file, "CELL_DATA %i\n", (iu-il+1)*(ju-jl+1)*(ku-kl+1));
    } else {
        fprintf(file, "CELL_DATA %i\n", (iu-il+1)*(ju-jl+1)*(ku-kl+1));
    }
    fprintf(file, "SCALARS pressure float 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");
    if (isBinaryVtk) {
        if (isMpiMode) {
            for (int k = kl; k <= ku; k++) {
                for (int j = jl; j <= ju; j++) {
                    for (int i = il; i <= iu; i++) {
                        if (isFluid(Flag[IDXFLAG_MPI(i,j,k)])) {
                            cellData[IDXCELL(i,j,k)] = P[IDXP_MPI(i,j,k)];
                        } else {
                            cellData[IDXCELL(i,j,k)] = 0;
                        }
                    }
                }
            }
        } else {
            #pragma omp parallel for shared(P, Flag) default(none)
            for (int k = kl; k <= ku; k++) {
                for (int j = jl; j <= ju; j++) {
                    for (int i = il; i <= iu; i++) {
                        if (isFluid(Flag[IDXFLAG_NORMAL(i,j,k)])) {
                            cellData[IDXCELL(i,j,k)] = P[IDXP_NORMAL(i,j,k)];
                        } else {
                            cellData[IDXCELL(i,j,k)] = 0;
                        }
                    }
                }
            }
        }
        swapEndianness(cellData, (iu-il+1)*(ju-jl+1)*(ku-kl+1));
        fwrite(cellData, sizeof(float), (iu-il+1)*(ju-jl+1)*(ku-kl+1), file);
    } else {
        if (isMpiMode) {
            for (int k = kl; k <= ku; k++) {
                for (int j = jl; j <= ju; j++) {
                    for (int i = il; i <= iu; i++) {
                        if (isFluid(Flag[IDXFLAG_MPI(i,j,k)])) {
                            fprintf(file, "%f\n", P[IDXP_MPI(i,j,k)]);
                        } else {
                            fprintf(file, "0\n");
                        }
                    }
                }
            }
        } else {
            for (int k = kl; k <= ku; k++) {
                for (int j = jl; j <= ju; j++) {
                    for (int i = il; i <= iu; i++) {
                        if (isFluid(Flag[IDXFLAG_NORMAL(i,j,k)])) {
                            fprintf(file, "%f\n", P[IDXP_NORMAL(i,j,k)]);
                        } else {
                            fprintf(file, "0\n");
                        }
                    }
                }
            }
        }
    }

    fprintf(file, "SCALARS temperature float 1\n");
    fprintf(file, "LOOKUP_TABLE default\n");
    if (isBinaryVtk) {
        if (isMpiMode) {
            for (int k = kl; k <= ku; k++) {
                for (int j = jl; j <= ju; j++) {
                    for (int i = il; i <= iu; i++) {
                        if (isFluid(Flag[IDXFLAG_MPI(i,j,k)])) {
                            cellData[IDXCELL(i,j,k)] = T[IDXT_MPI(i,j,k)];
                        } else {
                            cellData[IDXCELL(i,j,k)] = 0;
                        }
                    }
                }
            }
        } else {
            #pragma omp parallel for shared(T, Flag) default(none)
            for (int k = kl; k <= ku; k++) {
                for (int j = jl; j <= ju; j++) {
                    for (int i = il; i <= iu; i++) {
                        if (isFluid(Flag[IDXFLAG_NORMAL(i,j,k)])) {
                            cellData[IDXCELL(i,j,k)] = T[IDXT_NORMAL(i,j,k)];
                        } else {
                            cellData[IDXCELL(i,j,k)] = 0;
                        }
                    }
                }
            }
        }
        swapEndianness(cellData, (iu-il+1)*(ju-jl+1)*(ku-kl+1));
        fwrite(cellData, sizeof(float), (iu-il+1)*(ju-jl+1)*(ku-kl+1), file);
    } else {
        if (isMpiMode) {
            for (int k = kl; k <= ku; k++) {
                for (int j = jl; j <= ju; j++) {
                    for (int i = il; i <= iu; i++) {
                        if (isFluid(Flag[IDXFLAG_MPI(i,j,k)])) {
                            fprintf(file, "%f\n", T[IDXT_MPI(i,j,k)]);
                        } else {
                            fprintf(file, "0\n");
                        }
                    }
                }
            }
        } else {
            for (int k = kl; k <= ku; k++) {
                for (int j = jl; j <= ju; j++) {
                    for (int i = il; i <= iu; i++) {
                        if (isFluid(Flag[IDXFLAG_NORMAL(i,j,k)])) {
                            fprintf(file, "%f\n", T[IDXT_NORMAL(i,j,k)]);
                        } else {
                            fprintf(file, "0\n");
                        }
                    }
                }
            }
        }
    }

    if (isBinaryVtk) {
        fprintf(file, "SCALARS geometry unsigned_char 1\n");
        fprintf(file, "LOOKUP_TABLE default\n");
        if (isMpiMode) {
            for (int k = kl; k <= ku; k++) {
                for (int j = jl; j <= ju; j++) {
                    for (int i = il; i <= iu; i++) {
                        cellDataUint[IDXCELL(i,j,k)] = (Flag[IDXFLAG_MPI(i,j,k)] & 0x1) == 1;
                    }
                }
            }
        } else {
            #pragma omp parallel for shared(Flag) default(none)
            for (int k = kl; k <= ku; k++) {
                for (int j = jl; j <= ju; j++) {
                    for (int i = il; i <= iu; i++) {
                        cellDataUint[IDXCELL(i,j,k)] = (Flag[IDXFLAG_NORMAL(i,j,k)] & 0x1) == 1;
                    }
                }
            }
        }
        fwrite(cellDataUint, sizeof(uint8_t), (iu-il+1)*(ju-jl+1)*(ku-kl+1), file);
    } else {
        fprintf(file, "SCALARS geometry bit 1\n");
        fprintf(file, "LOOKUP_TABLE default\n");
        if (isMpiMode) {
            for (int k = kl; k <= ku; k++) {
                for (int j = jl; j <= ju; j++) {
                    for (int i = il; i <= iu; i++) {
                        fprintf(file, "%u\n", (Flag[IDXFLAG_MPI(i, j, k)] & 0x1) == 1);
                    }
                }
            }
        } else {
            for (int k = kl; k <= ku; k++) {
                for (int j = jl; j <= ju; j++) {
                    for (int i = il; i <= iu; i++) {
                        fprintf(file, "%u\n", (Flag[IDXFLAG_NORMAL(i, j, k)] & 0x1) == 1);
                    }
                }
            }
        }
    }
}
