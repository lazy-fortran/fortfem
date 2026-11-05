# Parametric P2/P1 BEM cost and normal-field audit

`fortfem_parametric_bem_scaling` compares full dense assembly with the five-sector
block-circulant assembler on an analytic corrugated torus. Its major radius is
5.5, minor radius 0.55 and helical ripple amplitude 0.1, with

    R = 5.5 + 0.55*cos(theta) + 0.1*cos(theta-5*phi),
    Z = 0.55*sin(theta) + 0.1*sin(theta-5*phi).

This synthetic five-period surface has no machine-data input. It is a simple
aspect-ratio and periodicity study, rather than a W7-X geometry surrogate with
validated machine accuracy. FortSym generates the position and angle derivatives.
The exterior harmonic oracle is the sum of five unit monopoles located inside
the tube at R=5.5, Z=0 and phi=2*pi*k/5. Its boundary normal derivative is known
independently of the BEM solve.

Run headlessly, choosing the thread count for the available resource allocation:

```sh
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 fpm run --profile release \
  --flag '-fopenmp -finline-matmul-limit=0' \
  --target fortfem_parametric_bem_scaling -- 4 40 circulant 1e-4 24
```

Arguments are the poloidal count, toroidal count, `dense`/`circulant`/`projection`,
pair relative tolerance, and audit quadrature degree. Toroidal counts must be a
multiple of five. There are twice the product of the two mesh counts panels.
The pair absolute tolerance is 1e-11 and maximum degree 64. A nonconvergent
assembly or solve stops explicitly. Refine the pair tolerance and audit degree
independently before interpreting a mesh ladder. A strong-scaling comparison
uses identical arguments with 1, 2, 4 and 8 OpenMP threads.

`projection` reports the quadrature-estimated best discontinuous P1 approximation
of the exact normal flux in the relative surface L2 norm. This isolates the
trace-space limitation and estimates a necessary mesh size. The solved error
also includes Dirichlet interpolation, quadrature and BEM errors. Neither norm
is a total three-dimensional field-error certificate. The audit integration is
independent of the quadrature accepted by the assembler.

Output reports assembly wall time, factorization plus one flux solve, both
normal-flux errors, stored operator bytes, pair work and degree, and Linux peak
resident memory when `/proc/self/status` is available. The solve timing includes
the right-hand-side application and factors afresh. It does not measure building
or applying a cached full DtN matrix. Process peak memory includes quadrature
samples and solver temporaries in addition to the reported operator storage.

On this genus-one mesh there are 3*N flux and 2*N potential DOFs. Dense V, K,
K-adjoint, W and mass arrays alone require 248*N^2 bytes, whereas the four stored
circulant block rows require 200*N^2/5 bytes. Dense expansion, Fourier matrices,
quadrature banks and factorization require additional memory. Compare the
measured costs and extrapolation limits with the chosen practical threshold
(64 GB or one hour for one assembly on an sCluster node) before selecting a
compressed operator.
