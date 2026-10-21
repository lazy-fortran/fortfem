#include "fortfem_element_cuda.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <type_traits>

static void check(cudaError_t result) {
    if (result != cudaSuccess) {
        std::fprintf(stderr,"CUDA: %s\n",cudaGetErrorString(result));
        std::exit(1);
    }
}
template<class S, class C> bool run() {
    constexpr int n=257;
    S *jac, *diffusion, *reaction;
    C *blocks;
    int *status;
    check(cudaMallocManaged(&jac,4*n*sizeof(S)));
    check(cudaMallocManaged(&diffusion,n*sizeof(S)));
    check(cudaMallocManaged(&reaction,n*sizeof(S)));
    check(cudaMallocManaged(&blocks,9*n*sizeof(C)));
    check(cudaMallocManaged(&status,n*sizeof(int)));
    for (int e=0; e<n; ++e) {
        jac[e]=S(1.0+double(e)/n); jac[e+n]=S(0.125);
        jac[e+2*n]=S(-0.25); jac[e+3*n]=S(0.75);
        if (e%2) { jac[e+2*n]=-jac[e+2*n]; jac[e+3*n]=-jac[e+3*n]; }
        diffusion[e]=S(0.5); reaction[e]=S(0.25);
    }
    // Singular final cell must be rejected before evaluating reciprocals.
    for (int j=0; j<4; ++j) jac[n-1+n*j]=0;
    if constexpr (std::is_same_v<S,double>) {
        check(fortfem_p1_batch_fp64(n,jac,diffusion,reaction,blocks,status,0));
    } else if constexpr (std::is_same_v<C,float>) {
        check(fortfem_p1_batch_fp32(n,jac,diffusion,reaction,blocks,status,0));
    } else {
        check(fortfem_p1_batch_mixed(n,jac,diffusion,reaction,blocks,status,0));
    }
    check(cudaDeviceSynchronize());
    bool passed=status[n-1]==2;
    for (int e=0; e<n-1; ++e) {
        double a=jac[e], b=jac[e+n], c=jac[e+2*n], d=jac[e+3*n];
        double det=a*d-b*c, area=std::abs(det)/2;
        double gx[3]={(b-d)/det,d/det,-b/det};
        double gy[3]={(c-a)/det,-c/det,a/det};
        double kappa=std::pow(std::abs(a)+std::abs(b)+std::abs(c)+std::abs(d),2)
            /std::abs(det);
        double u=std::numeric_limits<C>::epsilon()/2;
        double gamma=96*u/(1-96*u), reference[9], scale=0;
        for (int j=0; j<3; ++j) for (int i=0; i<3; ++i) {
            reference[i+3*j]=double(diffusion[e])*area*(gx[i]*gx[j]+gy[i]*gy[j])
                +double(reaction[e])*area*(i==j ? 2.0 : 1.0)/12;
            scale=std::fmax(scale,std::abs(reference[i+3*j]));
        }
        passed=passed && status[e]==0;
        for (int j=0; j<9; ++j)
            passed=passed && std::abs(double(blocks[e+n*j])-reference[j])
                <=gamma*kappa*kappa*scale;
    }
    for (int j=0; j<9; ++j) passed=passed && blocks[n-1+n*j]==0;
    check(cudaFree(jac)); check(cudaFree(diffusion)); check(cudaFree(reaction));
    check(cudaFree(blocks)); check(cudaFree(status));
    return passed;
}
int main() {
    int devices=0;
    if (cudaGetDeviceCount(&devices)!=cudaSuccess || devices==0) return 77;
    bool passed=run<double,double>() && run<float,float>() && run<float,double>();
    check(fortfem_p1_batch_fp64(0,nullptr,nullptr,nullptr,nullptr,nullptr,0));
    std::puts(passed ? "PASS CUDA cell roundoff oracles" : "FAIL CUDA cells");
    return passed ? 0 : 1;
}
