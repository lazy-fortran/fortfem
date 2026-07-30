# Adaptive Laplace BEM on a prolate spheroid

This example solves the exterior electrostatic Dirichlet problem with constant
P0 surface elements. It estimates the fine-grid Galerkin defect, marks the
smallest set of panels carrying half of its squared norm, performs conforming
red/green refinement, and projects new vertices back to the spheroid.

For semi-axes \(a>b\), the capacity in the kernel convention
\(G(x,y)=1/(4\pi|x-y|)\) is

\[
C=4\pi\frac{\sqrt{a^2-b^2}}{\operatorname{acosh}(a/b)}.
\]

The generated convergence plot compares the relative error against this
analytical value with the residual indicator. Generated images remain build
artifacts and are not checked into the repository.

Provenance:

- [Feischl et al., *Convergence of adaptive BEM and adaptive FEM-BEM
  coupling for estimators without h-weighting
  factor*](https://arxiv.org/abs/1405.5306) documents the
  solve-estimate-mark-refine loop, two-level indicators, and Dörfler marking.
- [Kraniotis and Leontaris, *Closed form solution for the surface area, the
  capacitance and the demagnetizing factors of the
  ellipsoid*](https://arxiv.org/abs/1306.0509) gives the independent prolate
  spheroid capacity formula. The factor \(4\pi\) converts their
  sphere-capacity convention to FortFEM's \(1/(4\pi r)\) kernel.

Run with:

```sh
fpm run --example adaptive_bem_prolate
```
