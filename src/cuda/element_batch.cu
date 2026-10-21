#include "fortfem_element_cuda.h"
#include "fortfem_element_precision_products.cuh"
#include <limits>
#include <type_traits>

template<class S, class C>
__global__ void p1_batch(int n, const S *jac, const S *diffusion,
    const S *reaction, C *blocks, int *status, C epsilon) {
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= n) return;
    C a = jac[e], b = jac[e+n], c = jac[e+2*n], d = jac[e+3*n];
    C k = diffusion[e], m = reaction[e], det, scale, local[9] = {};
    int code = 1;
    if (isfinite(a) && isfinite(b) && isfinite(c) && isfinite(d) &&
        isfinite(k) && isfinite(m)) {
        if constexpr (std::is_same_v<C, float>)
            generated_triangle_measure_fp32(a,b,c,d,&det,&scale);
        else
            generated_triangle_measure_fp64(a,b,c,d,&det,&scale);
        code = 2;
        if (isfinite(scale) && fabs(det) > epsilon*scale) {
            if constexpr (std::is_same_v<C, float>)
                generated_p1_cell_fp32(a,b,c,d,k,m,local,local+1,local+2,
                    local+3,local+4,local+5,local+6,local+7,local+8);
            else
                generated_p1_cell_fp64(a,b,c,d,k,m,local,local+1,local+2,
                    local+3,local+4,local+5,local+6,local+7,local+8);
            code = 0;
            for (int j=0; j<9; ++j) if (!isfinite(local[j])) code = 3;
        }
    }
    status[e] = code;
    for (int j=0; j<9; ++j) blocks[e+n*j] = code == 0 ? local[j] : C(0);
}

template<class S, class C>
cudaError_t launch(int n, const S *jac, const S *diffusion, const S *reaction,
    C *blocks, int *status, cudaStream_t stream) {
    if (n < 0) return cudaErrorInvalidValue;
    if (n == 0) return cudaSuccess;
    p1_batch<S,C><<<(n+255)/256,256,0,stream>>>(
        n,jac,diffusion,reaction,blocks,status,std::numeric_limits<C>::epsilon());
    return cudaGetLastError();
}
extern "C" cudaError_t fortfem_p1_batch_fp32(int n, const float *jac,
    const float *diffusion, const float *reaction, float *blocks,
    int *status, cudaStream_t stream) {
    return launch(n,jac,diffusion,reaction,blocks,status,stream);
}
extern "C" cudaError_t fortfem_p1_batch_fp64(int n, const double *jac,
    const double *diffusion, const double *reaction, double *blocks,
    int *status, cudaStream_t stream) {
    return launch(n,jac,diffusion,reaction,blocks,status,stream);
}
extern "C" cudaError_t fortfem_p1_batch_mixed(int n, const float *jac,
    const float *diffusion, const float *reaction, double *blocks,
    int *status, cudaStream_t stream) {
    return launch(n,jac,diffusion,reaction,blocks,status,stream);
}
