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
 
kernel void calculateFghOpencl(
        Real Re, Real GX, Real GY, Real GZ, Real alpha, Real beta,
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        global Real *U, global Real *V, global Real *W, global Real *T, global Real *F, global Real *G, global Real *H,
        global FlagType *Flag) {
    int i = get_global_id(2) + 1;
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    Real d2u_dx2,d2u_dy2,d2u_dz2,
            d2v_dx2,d2v_dy2,d2v_dz2,
            d2w_dx2,d2w_dy2,d2w_dz2;

    Real du2_dx,duv_dy,duw_dz,
            duv_dx,dv2_dy,dvw_dz,
            duw_dx,dvw_dy,dw2_dz;

    Real Dx = 1/dx, Dy = 1/dy, Dz = 1/dz;

    if (i <= imax - 1 && j <= jmax && k <= kmax) {
        if(isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i+1,j,k)])){
            d2u_dx2 = (U[IDXU(i+1,j,k)] - 2*U[IDXU(i,j,k)] + U[IDXU(i-1,j,k)])/(dx*dx);
            d2u_dy2 = (U[IDXU(i,j+1,k)] - 2*U[IDXU(i,j,k)] + U[IDXU(i,j-1,k)])/(dy*dy);
            d2u_dz2 = (U[IDXU(i,j,k+1)] - 2*U[IDXU(i,j,k)] + U[IDXU(i,j,k-1)])/(dz*dz);

            du2_dx = (Real)(0.25)*Dx*(
                    (U[IDXU(i,j,k)]+U[IDXU(i+1,j,k)])*(U[IDXU(i,j,k)]+U[IDXU(i+1,j,k)]) -
                    (U[IDXU(i-1,j,k)]+U[IDXU(i,j,k)])*(U[IDXU(i-1,j,k)]+U[IDXU(i,j,k)]) +
                    alpha*(
                            (fabs(U[IDXU(i,j,k)]+U[IDXU(i+1,j,k)])*(U[IDXU(i,j,k)]-U[IDXU(i+1,j,k)]))-
                            (fabs(U[IDXU(i-1,j,k)]+U[IDXU(i,j,k)])*(U[IDXU(i-1,j,k)]-U[IDXU(i,j,k)]))
                    )
            );

            duv_dy = (Real)(0.25)*Dy*(
                    (V[IDXV(i,j,k)]+V[IDXV(i+1,j,k)])*(U[IDXU(i,j,k)]+U[IDXU(i,j+1,k)]) -
                    (V[IDXV(i,j-1,k)]+V[IDXV(i+1,j-1,k)])*(U[IDXU(i,j-1,k)]+U[IDXU(i,j,k)]) +
                    alpha*(
                            (fabs(V[IDXV(i,j,k)]+V[IDXV(i+1,j,k)])*(U[IDXU(i,j,k)]-U[IDXU(i,j+1,k)]))-
                            (fabs(V[IDXV(i,j-1,k)]+V[IDXV(i+1,j-1,k)])*(U[IDXU(i,j-1,k)]-U[IDXU(i,j,k)]))
                    )
            );

            duw_dz = (Real)(0.25)*Dz*(
                    (W[IDXW(i,j,k)]+W[IDXW(i+1,j,k)])*(U[IDXU(i,j,k)]+U[IDXU(i,j,k+1)]) -
                    (W[IDXW(i,j,k-1)]+W[IDXW(i+1,j,k-1)])*(U[IDXU(i,j,k-1)]+U[IDXU(i,j,k)]) +
                    alpha*(
                            (fabs(W[IDXW(i,j,k)]+W[IDXW(i+1,j,k)])*(U[IDXU(i,j,k)]-U[IDXU(i,j,k+1)]))-
                            (fabs(W[IDXW(i,j,k-1)]+W[IDXW(i+1,j,k-1)])*(U[IDXU(i,j,k-1)]-U[IDXU(i,j,k)]))
                    )
            );

            F[IDXF(i,j,k)] = U[IDXU(i,j,k)] + dt * (
                    (1/Re)*(d2u_dx2+d2u_dy2+d2u_dz2)-
                    du2_dx-duv_dy-duw_dz+
                    GX-(beta/2)*(T[IDXT(i,j,k)]+T[IDXT(i+1,j,k)])*GX
            );
        } else if (B_L(Flag[IDXFLAG(i,j,k)])) {
                F[IDXF(i-1,j,k)] = U[IDXU(i-1,j,k)];
        } else if (B_R(Flag[IDXFLAG(i,j,k)])) {
                F[IDXF(i,j,k)] = U[IDXU(i,j,k)];
        }
    }

    if (i <= imax && j <= jmax - 1 && k <= kmax){
        if(isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i,j+1,k)])){
            d2v_dx2 = (V[IDXV(i+1,j,k)] - 2*V[IDXV(i,j,k)] + V[IDXV(i-1,j,k)])/(dx*dx);
            d2v_dy2 = (V[IDXV(i,j+1,k)] - 2*V[IDXV(i,j,k)] + V[IDXV(i,j-1,k)])/(dy*dy);
            d2v_dz2 = (V[IDXV(i,j,k+1)] - 2*V[IDXV(i,j,k)] + V[IDXV(i,j,k-1)])/(dz*dz);

            duv_dx = (Real)(0.25)*Dx*(
                    (U[IDXU(i,j,k)]+U[IDXU(i,j+1,k)])*(V[IDXV(i,j,k)]+V[IDXV(i+1,j,k)]) -
                    (U[IDXU(i-1,j,k)]+U[IDXU(i-1,j+1,k)])*(V[IDXV(i-1,j,k)]+V[IDXV(i,j,k)]) +
                    alpha*(
                            (fabs(U[IDXU(i,j,k)]+U[IDXU(i,j+1,k)])*(V[IDXV(i,j,k)]-V[IDXV(i+1,j,k)]))-
                            (fabs(U[IDXU(i-1,j,k)]+U[IDXU(i-1,j+1,k)])*(V[IDXV(i-1,j,k)]-V[IDXV(i,j,k)]))
                    )
            );

            dv2_dy = (Real)(0.25)*Dy*(
                    (V[IDXV(i,j,k)]+V[IDXV(i,j+1,k)])*(V[IDXV(i,j,k)]+V[IDXV(i,j+1,k)]) -
                    (V[IDXV(i,j-1,k)]+V[IDXV(i,j,k)])*(V[IDXV(i,j-1,k)]+V[IDXV(i,j,k)]) +
                    alpha*(
                            (fabs(V[IDXV(i,j,k)]+V[IDXV(i,j+1,k)])*(V[IDXV(i,j,k)]-V[IDXV(i,j+1,k)]))-
                            (fabs(V[IDXV(i,j-1,k)]+V[IDXV(i,j,k)])*(V[IDXV(i,j-1,k)]-V[IDXV(i,j,k)]))
                    )
            );

            dvw_dz = (Real)(0.25)*Dz*(
                    (W[IDXW(i,j,k)]+W[IDXW(i,j+1,k)])*(V[IDXV(i,j,k)]+V[IDXV(i,j,k+1)]) -
                    (W[IDXW(i,j,k-1)]+W[IDXW(i,j+1,k-1)])*(V[IDXV(i,j,k-1)]+V[IDXV(i,j,k)]) +
                    alpha*(
                            (fabs(W[IDXW(i,j,k)]+W[IDXW(i,j+1,k)])*(V[IDXV(i,j,k)]-V[IDXV(i,j,k+1)]))-
                            (fabs(W[IDXW(i,j,k-1)]+W[IDXW(i,j+1,k-1)])*(V[IDXV(i,j,k-1)]-V[IDXV(i,j,k)]))
                    )
            );

            G[IDXG(i,j,k)] = V[IDXV(i,j,k)] + dt * (
                    (1/Re)*(d2v_dx2+d2v_dy2+d2v_dz2)-
                    duv_dx-dv2_dy-dvw_dz+
                    GY-(beta/2)*(T[IDXT(i,j,k)]+T[IDXT(i,j+1,k)])*GY
            );
        } else if (B_D(Flag[IDXFLAG(i,j,k)])) {
            G[IDXG(i,j-1,k)] = V[IDXV(i,j-1,k)];
        } else if (B_U(Flag[IDXFLAG(i,j,k)])) {
            G[IDXG(i,j,k)] = V[IDXV(i,j,k)];
        }
    }

    if (i <= imax && j <= jmax && k <= kmax - 1) {
        if (isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i,j,k+1)])){
            d2w_dx2 = (W[IDXW(i+1,j,k)] - 2*W[IDXW(i,j,k)] + W[IDXW(i-1,j,k)])/(dx*dx);
            d2w_dy2 = (W[IDXW(i,j+1,k)] - 2*W[IDXW(i,j,k)] + W[IDXW(i,j-1,k)])/(dy*dy);
            d2w_dz2 = (W[IDXW(i,j,k+1)] - 2*W[IDXW(i,j,k)] + W[IDXW(i,j,k-1)])/(dz*dz);

            duw_dx = (Real)(0.25)*Dx*(
                    (U[IDXU(i,j,k)]+U[IDXU(i,j,k+1)])*(W[IDXW(i,j,k)]+W[IDXW(i+1,j,k)]) -
                    (U[IDXU(i-1,j,k)]+U[IDXU(i-1,j,k+1)])*(W[IDXW(i-1,j,k)]+W[IDXW(i,j,k)]) +
                    alpha*(
                            (fabs(U[IDXU(i,j,k)]+U[IDXU(i,j,k+1)])*(W[IDXW(i,j,k)]-W[IDXW(i+1,j,k)]))-
                            (fabs(U[IDXU(i-1,j,k)]+U[IDXU(i-1,j,k+1)])*(W[IDXW(i-1,j,k)]-W[IDXW(i,j,k)]))
                    )
            );

            dvw_dy = (Real)(0.25)*Dy*(
                    (V[IDXV(i,j,k)]+V[IDXV(i,j,k+1)])*(W[IDXW(i,j,k)]+W[IDXW(i,j+1,k)]) -
                    (V[IDXV(i,j-1,k)]+V[IDXV(i,j-1,k+1)])*(W[IDXW(i,j-1,k)]+W[IDXW(i,j,k)]) +
                    alpha*(
                            (fabs(V[IDXV(i,j,k)]+V[IDXV(i,j,k+1)])*(W[IDXW(i,j,k)]-W[IDXW(i,j+1,k)]))-
                            (fabs(V[IDXV(i,j-1,k)]+V[IDXV(i,j-1,k+1)])*(W[IDXW(i,j-1,k)]-W[IDXW(i,j,k)]))
                    )
            );

            dw2_dz = (Real)(0.25)*Dz*(
                    (W[IDXW(i,j,k)]+W[IDXW(i,j,k+1)])*(W[IDXW(i,j,k)]+W[IDXW(i,j,k+1)]) -
                    (W[IDXW(i,j,k-1)]+W[IDXW(i,j,k)])*(W[IDXW(i,j,k-1)]+W[IDXW(i,j,k)]) +
                    alpha*(
                            (fabs(W[IDXW(i,j,k)]+W[IDXW(i,j,k+1)])*(W[IDXW(i,j,k)]-W[IDXW(i,j,k+1)]))-
                            (fabs(W[IDXW(i,j,k-1)]+W[IDXW(i,j,k)])*(W[IDXW(i,j,k-1)]-W[IDXW(i,j,k)]))
                    )
            );

            H[IDXH(i,j,k)] = W[IDXW(i,j,k)] +  dt * (
                    (1/Re)*(d2w_dx2+d2w_dy2+d2w_dz2)-
                    duw_dx-dvw_dy-dw2_dz+
                    GZ-(beta/2)*(T[IDXT(i,j,k)]+T[IDXT(i,j,k+1)])*GZ
            );
        } else if (B_B(Flag[IDXFLAG(i,j,k)])) {
            H[IDXH(i,j,k-1)] = W[IDXW(i,j,k-1)];
        } else if (B_F(Flag[IDXFLAG(i,j,k)])) {
            H[IDXH(i,j,k)] = W[IDXW(i,j,k)];
        }
    }
}

kernel void setFBoundariesOpencl(int imax, int jmax, int kmax, global Real *U, global Real *F) {
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    // Set the boundary values for F on the y-z-planes.
    if (j <= jmax && k <= kmax){
        F[IDXF(0,j,k)] = U[IDXU(0,j,k)];
        F[IDXF(imax,j,k)] = U[IDXU(imax,j,k)];
    }
}

kernel void setGBoundariesOpencl(int imax, int jmax, int kmax, global Real *V, global Real *G) {
    int i = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    // Set the boundary values for G on the x-z-planes.
    if (i <= imax && k <= kmax){
        G[IDXG(i,0,k)] = V[IDXV(i,0,k)];
        G[IDXG(i,jmax,k)] = V[IDXV(i,jmax,k)];
    }
}

kernel void setHBoundariesOpencl(int imax, int jmax, int kmax, global Real *W, global Real *H) {

    int i = get_global_id(1) + 1;
    int j = get_global_id(0) + 1;

    // Set the boundary values for G on the x-z-planes.
    if (i <= imax && j <= jmax){
        H[IDXH(i,j,0)] = W[IDXW(i,j,0)];
        H[IDXH(i,j,kmax)] = W[IDXW(i,j,kmax)];
    }
}

kernel void calculateRsOpencl(
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        global Real *F, global Real *G, global Real *H, global Real *RS) {
    int i = get_global_id(2) + 1;
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    if (i <= imax && j <= jmax && k <= kmax){
        RS[IDXRS(i, j, k)] = ((F[IDXF(i, j, k)] - F[IDXF(i - 1, j, k)]) / dx +
                (G[IDXG(i, j, k)] - G[IDXG(i, j - 1, k)]) / dy +
                (H[IDXH(i, j, k)] - H[IDXH(i, j, k - 1)]) / dz) / dt;
    }
}

/**
 * Reference: Based on kernel 4 from https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
 * @param input The array of input values (of size 'sizeOfInput').
 * @param output The output array (of size iceil(numberOfBlocksI, blockSize1D*2)).
 * @param sdata The local memory to use (of byte size blockSize1D * sizeof(Real)).
 * @param sizeOfInput The number of input values.
 */
kernel void calculateMaximum(
        global Real *input, global Real *output, local Real *sdata, int sizeOfInput) {
    unsigned int threadID = get_local_id(0);
    unsigned int i = get_group_id(0) * get_local_size(0) * 2 + get_local_id(0);

    // Copy the data to the shared memory and do the first reduction step.
    if (i + get_local_size(0) < sizeOfInput){
        sdata[threadID] = fmax(fabs(input[i]), fabs(input[i + get_local_size(0)]));
    } else if (i < sizeOfInput){
        sdata[threadID] = input[i];
    } else{
        sdata[threadID] = 0;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // Do the reduction in the shared memory.
    for (unsigned int stride = get_local_size(0) / 2; stride > 0; stride >>= 1) {
        if (threadID < stride) {
            sdata[threadID] = fmax(fabs(sdata[threadID]), fabs(sdata[threadID + stride]));
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // Write the result for this block to global memory.
    if (threadID == 0) {
        output[get_group_id(0)] = sdata[0];
    }
}

kernel void calculateUvwOpencl(
        Real dt, Real dx, Real dy, Real dz, int imax, int jmax, int kmax,
        global Real *U, global Real *V, global Real *W, global Real *F, global Real *G, global Real *H, global Real *P,
        global FlagType *Flag) {
    int i = get_global_id(2) + 1;
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    if (i <= imax - 1 && j <= jmax && k <= kmax){
        if(isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i+1,j,k)])){
            U[IDXU(i, j, k)] = F[IDXF(i, j, k)] - dt / dx * (P[IDXP(i + 1, j, k)] - P[IDXP(i, j, k)]);
        }
    }

    if (i <= imax && j <= jmax - 1 && k <= kmax){
        if(isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i,j+1,k)])){
            V[IDXV(i, j, k)] = G[IDXG(i, j, k)] - dt / dy * (P[IDXP(i, j + 1, k)] - P[IDXP(i, j, k)]);
        }
    }

    if (i <= imax && j <= jmax && k <= kmax - 1){
        if(isFluid(Flag[IDXFLAG(i,j,k)]) && isFluid(Flag[IDXFLAG(i,j,k+1)])){
            W[IDXW(i, j, k)] = H[IDXH(i, j, k)] - dt / dz * (P[IDXP(i, j, k + 1)] - P[IDXP(i, j, k)]);
        }
    }
}

kernel void calculateTemperatureOpencl(
        Real Re, Real Pr, Real alpha,
        Real dt, Real dx, Real dy, Real dz,
        int imax, int jmax, int kmax,
        global Real *U, global Real *V, global Real *W, global Real *T, global Real *T_temp, global FlagType *Flag) {
    int i = get_global_id(2) + 1;
    int j = get_global_id(1) + 1;
    int k = get_global_id(0) + 1;

    Real duT_dx, dvT_dy, dwT_dz, d2T_dx2, d2T_dy2, d2T_dz2;

    if (i <= imax && j <= jmax && k <= kmax){
        if(isFluid(Flag[IDXFLAG(i,j,k)])){
            duT_dx = 1 / dx * (
                    U[IDXU(i, j, k)] * ((T_temp[IDXT(i, j, k)] + T_temp[IDXT(i + 1, j, k)]) / 2) -
                    U[IDXU(i - 1, j, k)] * ((T_temp[IDXT(i - 1, j, k)] + T_temp[IDXT(i, j, k)]) / 2) +
                    alpha * (
                            fabs(U[IDXU(i, j, k)])*((T_temp[IDXT(i, j, k)] - T_temp[IDXT(i + 1, j, k)]) / 2) -
                            fabs(U[IDXU(i - 1, j, k)])*((T_temp[IDXT(i - 1, j, k)] - T_temp[IDXT(i, j, k)]) / 2)
                    )
            );

            dvT_dy = 1 / dy * (
                    V[IDXV(i, j, k)] * ((T_temp[IDXT(i, j, k)] + T_temp[IDXT(i, j + 1, k)]) / 2) -
                    V[IDXV(i, j - 1, k)] * ((T_temp[IDXT(i, j - 1, k)] + T_temp[IDXT(i, j, k)]) / 2) +
                    alpha * (
                            fabs(V[IDXV(i, j, k)])*((T_temp[IDXT(i, j, k)] - T_temp[IDXT(i, j + 1, k)]) / 2) -
                            fabs(V[IDXV(i, j - 1, k)])*((T_temp[IDXT(i, j - 1, k)] - T_temp[IDXT(i, j, k)]) / 2)
                    )
            );

            dwT_dz = 1 / dz * (
                    W[IDXW(i, j, k)] * ((T_temp[IDXT(i, j, k)] + T_temp[IDXT(i, j, k + 1)]) / 2) -
                    W[IDXW(i, j, k - 1)] * ((T_temp[IDXT(i, j, k - 1)] + T_temp[IDXT(i, j, k)]) / 2) +
                    alpha * (
                            fabs(W[IDXW(i, j, k)])*((T_temp[IDXT(i, j, k)] - T_temp[IDXT(i, j, k + 1)]) / 2) -
                            fabs(W[IDXW(i, j, k - 1)])*((T_temp[IDXT(i, j, k - 1)] - T_temp[IDXT(i, j, k)]) / 2)
                    )
            );

            d2T_dx2 =
                    (T_temp[IDXT(i + 1, j, k)] - 2 * T_temp[IDXT(i, j, k)] + T_temp[IDXT(i - 1, j, k)]) / (dx*dx);

            d2T_dy2 =
                    (T_temp[IDXT(i, j + 1, k)] - 2 * T_temp[IDXT(i, j, k)] + T_temp[IDXT(i, j - 1, k)]) / (dy*dy);

            d2T_dz2 =
                    (T_temp[IDXT(i, j, k + 1)] - 2 * T_temp[IDXT(i, j, k)] + T_temp[IDXT(i, j, k - 1)]) / (dz*dz);

            T[IDXT(i, j, k)] = T_temp[IDXT(i, j, k)] + dt * (
                    (1 / (Re*Pr))*(d2T_dx2 + d2T_dy2 + d2T_dz2) -
                    duT_dx -
                    dvT_dy -
                    dwT_dz
            );
        }
    }
}
