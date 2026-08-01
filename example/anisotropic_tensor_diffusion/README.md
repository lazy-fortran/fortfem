# Anisotropic tensor diffusion

This example solves a manufactured two-dimensional diffusion problem with a
strongly anisotropic constant tensor,

\[
-\nabla\!\cdot(K\nabla u)=f,
\qquad K=\operatorname{diag}(1000,1),
\qquad u=\sin(\pi x)\sin(\pi y).
\]

The P1 triangle assembly calls the public tensor diffusion contraction for
every element. The source is evaluated from the exact solution, and the
program checks the resulting numerical field, energy positivity, and
manufactured error before writing plots. The first plot is the computed
solution; the centerline and mesh views are secondary diagnostics.
An absolute-error map is included alongside the centerline and benchmark
values.
