# FreeFEM Poisson interoperability runner

This runner solves the shared unit-square manufactured problem

\[
-\Delta u=f,\qquad
u=\sin(\pi x)\sin(\pi y),\qquad
f=2\pi^2u,
\]

with homogeneous Dirichlet data and triangular P1 elements.

```sh
FreeFem++ poisson_2d.edp
```

It writes `poisson.csv` and exits nonzero if either analytical error does not
decrease. Generated records are CI artifacts rather than repository content.
