#pragma once
#include <cuda_runtime_api.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Device-resident SoA: jac[e+n*c], blocks[e+n*(row+3*column)].
 * Launches on the supplied stream; ownership and synchronization stay with caller.
 * status per cell: 0 success, 1 nonfinite input, 2 geometry, 3 overflow.
 */
cudaError_t fortfem_p1_batch_fp32(int n, const float *jac,
    const float *diffusion, const float *reaction, float *blocks,
    int *status, cudaStream_t stream);
cudaError_t fortfem_p1_batch_fp64(int n, const double *jac,
    const double *diffusion, const double *reaction, double *blocks,
    int *status, cudaStream_t stream);
cudaError_t fortfem_p1_batch_mixed(int n, const float *jac,
    const float *diffusion, const float *reaction, double *blocks,
    int *status, cudaStream_t stream);
#ifdef __cplusplus
}
#endif
