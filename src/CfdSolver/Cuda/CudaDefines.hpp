//
// Created by christoph on 10.07.19.
//

#ifndef CFD3D_CUDADEFINES_HPP
#define CFD3D_CUDADEFINES_HPP

__device__ inline bool isFluid(unsigned int flag) { return (flag >> 0) & 1; }
__device__ inline bool isNoSlip(unsigned int flag) { return (flag >> 1) & 1; }
__device__ inline bool isFreeSlip(unsigned int flag) { return (flag >> 2) & 1; }
__device__ inline bool isOutflow(unsigned int flag) { return (flag >> 3) & 1; }
__device__ inline bool isInflow(unsigned int flag) { return (flag >> 4) & 1; }
__device__ inline bool B_L(unsigned int flag) { return (flag >> 5) & 1; }
__device__ inline bool B_R(unsigned int flag) { return (flag >> 6) & 1; }
__device__ inline bool B_D(unsigned int flag) { return (flag >> 7) & 1; }
__device__ inline bool B_U(unsigned int flag) { return (flag >> 8) & 1; }
__device__ inline bool B_B(unsigned int flag) { return (flag >> 9) & 1; }
__device__ inline bool B_F(unsigned int flag) { return (flag >> 10) & 1; }
__device__ inline bool isHot(unsigned int flag) { return (flag >> 11) & 1; }
__device__ inline bool isCold(unsigned int flag) { return (flag >> 12) & 1; }
__device__ inline bool isCoupling(unsigned int flag) { return (flag >> 13) & 1; }

#ifndef TOSTRING
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#endif

void _checkCudaError(cudaError_t result, const char *expression, const char *locationText);
#define checkCudaError(expr) _checkCudaError(expr, #expr, __FILE__ ":" TOSTRING(__LINE__))

#endif //CFD3D_CUDADEFINES_HPP
