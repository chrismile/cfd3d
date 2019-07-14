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
 
kernel void setXYPlanesPressureBoundariesOpencl(
        int imax, int jmax, int kmax, global Real *P) {
    int i = get_global_id(1) + 1;
    int j = get_global_id(0) + 1;

    // Set the boundary values for the pressure on the x-y-planes.
    if (i <= imax && j <= jmax) {
        P[IDXP(i, j, 0)] = P[IDXP(i, j, 1)];
        P[IDXP(i, j, kmax + 1)] = P[IDXP(i, j, kmax)];
    }
}

kernel void setXZPlanesPressureBoundariesOpencl(
        int imax, int jmax, int kmax, global Real *P) {
    int i = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;
    // Set the boundary values for the pressure on the x-z-planes.
    if (i <= imax && k <= kmax) {
        P[IDXP(i, 0, k)] = P[IDXP(i, 1, k)];
        P[IDXP(i, jmax + 1, k)] = P[IDXP(i, jmax, k)];
    }
}

kernel void setYZPlanesPressureBoundariesOpencl(
        int imax, int jmax, int kmax, global Real *P) {
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;
    // Set the boundary values for the pressure on the y-z-planes.
    if (j <= jmax && k <= kmax) {
        P[IDXP(0, j, k)] = P[IDXP(1, j, k)];
        P[IDXP(imax + 1, j, k)] = P[IDXP(imax, j, k)];
    }
}

kernel void setBoundaryConditionsPressureInDomainOpencl(
        int imax, int jmax, int kmax, global Real *P, global FlagType *Flag) {
    int i = get_global_id(2) + 1;
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    if (i <= imax && j <= jmax && k <= kmax && !isFluid(Flag[IDXFLAG(i, j, k)])) {
        int numDirectFlag = 0;
        Real P_temp = (Real)(0);

        if (B_R(Flag[IDXFLAG(i, j, k)])) {
            P_temp += P[IDXP(i + 1, j, k)];
            numDirectFlag++;
        }

        if (B_L(Flag[IDXFLAG(i, j, k)])) {
            P_temp += P[IDXP(i - 1, j, k)];
            numDirectFlag++;
        }

        if (B_U(Flag[IDXFLAG(i, j, k)])) {
            P_temp += P[IDXP(i, j + 1, k)];
            numDirectFlag++;
        }

        if (B_D(Flag[IDXFLAG(i, j, k)])) {
            P_temp += P[IDXP(i, j - 1, k)];
            numDirectFlag++;
        }

        if (B_B(Flag[IDXFLAG(i, j, k)])) {
            P_temp += P[IDXP(i, j, k - 1)];
            numDirectFlag++;
        }

        if (B_F(Flag[IDXFLAG(i, j, k)])) {
            P_temp += P[IDXP(i, j, k + 1)];
            numDirectFlag++;
        }

        if (numDirectFlag == 0) {
            P[IDXP(i, j, k)] = 0;
        } else {
            P[IDXP(i, j, k)] = P_temp / (Real)(numDirectFlag);
        }
    }
}

kernel void copyPressureOpencl(int imax, int jmax, int kmax, global Real *P, global Real *P_temp) {
    int i = get_global_id(2);
    int j = get_global_id(1);
    int k = get_global_id(0);

    if (i <= imax + 1 && j <= jmax + 1 && k <= kmax + 1) {
        P_temp[IDXP(i, j, k)] = P[IDXP(i, j, k)];
    }
}

/**
 * Reference: Based on kernel 4 from https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
 * @param input The array of input values (of size 'sizeOfInput').
 * @param output The output array (of size iceil(numberOfBlocksI, blockSize1D*2)).
 * @param sdata The local memory to use (of byte size blockSize1D * sizeof(Real)).
 * @param sizeOfInput The number of input values.
 */
kernel void reduceSumOpenclKernelReal(global Real *input, global Real *output, local Real *sdata, int sizeOfInput) {
    unsigned int threadID = get_local_id(0);
    unsigned int i = get_group_id(0) * get_local_size(0) * 2 + get_local_id(0);

    // Copy the data to the shared memory and do the first reduction step.
    if (i + get_local_size(0) < sizeOfInput){
        sdata[threadID] = input[i] + input[i + get_local_size(0)];
    } else if (i < sizeOfInput){
        sdata[threadID] = input[i];
    } else{
        sdata[threadID] = 0;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // Do the reduction in the shared memory.
    for (unsigned int stride = get_local_size(0) / 2; stride > 0; stride >>= 1) {
        if (threadID < stride) {
            sdata[threadID] += sdata[threadID + stride];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // Write the result for this block to global memory.
    if (threadID == 0) {
        output[get_group_id(0)] = sdata[0];
    }
}

/**
 * Reference: Based on kernel 4 from https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
 * @param input The array of input values (of size 'sizeOfInput').
 * @param output The output array (of size iceil(numberOfBlocksI, blockSize1D*2)).
 * @param sdata The local memory to use (of byte size blockSize1D * sizeof(unsigned int)).
 * @param sizeOfInput The number of input values.
 */
kernel void reduceSumOpenclKernelUint(
        global unsigned int *input, global unsigned int *output, local unsigned int *sdata, int sizeOfInput) {
    unsigned int threadID = get_local_id(0);
    unsigned int i = get_group_id(0) * get_local_size(0) * 2 + get_local_id(0);

    // Copy the data to the shared memory and do the first reduction step.
    if (i + get_local_size(0) < sizeOfInput){
        sdata[threadID] = input[i] + input[i + get_local_size(0)];
    } else if (i < sizeOfInput){
        sdata[threadID] = input[i];
    } else{
        sdata[threadID] = 0;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // Do the reduction in the shared memory.
    for (unsigned int stride = get_local_size(0) / 2; stride > 0; stride >>= 1) {
        if (threadID < stride) {
            sdata[threadID] += sdata[threadID + stride];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // Write the result for this block to global memory.
    if (threadID == 0) {
        output[get_group_id(0)] = sdata[0];
    }
}

kernel void sorSolverIterationOpencl(
        Real omg, Real dx, Real dy, Real dz, Real coeff, int imax, int jmax, int kmax,
        global Real *P, global Real *P_temp, global Real *RS, global FlagType *Flag) {
    int i = get_global_id(2) + 1;
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    if (i <= imax && j <= jmax && k <= kmax) {
        if (isFluid(Flag[IDXFLAG(i,j,k)])) {
            P[IDXP(i, j, k)] = ((Real)(1.0) - omg) * P_temp[IDXP(i, j, k)] + coeff *
                    ((P_temp[IDXP(i + 1, j, k)] + P_temp[IDXP(i - 1, j, k)]) / (dx * dx)
                     + (P_temp[IDXP(i, j + 1, k)] + P_temp[IDXP(i, j - 1, k)]) / (dy * dy)
                     + (P_temp[IDXP(i, j, k + 1)] + P_temp[IDXP(i, j, k - 1)]) / (dz * dz)
                     - RS[IDXRS(i, j, k)]);
        }
    }
}

kernel void sorSolverComputeResidualArrayOpencl(
        Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        global Real *P, global Real *RS, global FlagType *Flag, global Real *residualArray,
        global unsigned int *numFluidCellsArray) {
    int i = get_global_id(2) + 1;
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    if (i <= imax && j <= jmax && k <= kmax){
        int arrayIndex1D = (i-1)*jmax*kmax + (j-1)*kmax + k-1;
        if (isFluid(Flag[IDXFLAG(i,j,k)])){
            residualArray[arrayIndex1D] = SQR(
                    (P[IDXP(i+1,j,k)] - (Real)(2.0)*P[IDXP(i,j,k)] + P[IDXP(i-1,j,k)])/(dx*dx)
                    + (P[IDXP(i,j+1,k)] - (Real)(2.0)*P[IDXP(i,j,k)] + P[IDXP(i,j-1,k)])/(dy*dy)
                    + (P[IDXP(i,j,k+1)] - (Real)(2.0)*P[IDXP(i,j,k)] + P[IDXP(i,j,k-1)])/(dz*dz)
                    - RS[IDXRS(i,j,k)]
            );
            numFluidCellsArray[arrayIndex1D] = 1;
        } else {
            residualArray[arrayIndex1D] = (Real)(0.0);
            numFluidCellsArray[arrayIndex1D] = 0;
        }
    }
}
