# Adaptive outgoing Helmholtz BEM on a sphere

This example applies the solve-estimate-mark-refine loop to the complex P0
single-layer equation for a unit sphere. New vertices are projected back to
the analytical surface after conforming red/green refinement.

For constant Dirichlet data and the outgoing convention used by FortFEM, the
radial exterior solution is

\[
u(r)=\frac{\exp(ik(r-1))}{r}.
\]

The field at \(r=2\) is therefore an independent analytical oracle. The
fine-grid Galerkin residual supplies local indicators, and Dörfler marking
selects panels carrying half of their squared norm.

Provenance:

- [Feischl et al., *Convergence of adaptive BEM and adaptive FEM-BEM
  coupling for estimators without h-weighting
  factor*](https://arxiv.org/abs/1405.5306) gives the adaptive loop and
  two-level-estimator setting.
- [NIST DLMF §10.47](https://dlmf.nist.gov/10.47) defines the outgoing
  spherical Hankel functions; its order-zero case reduces to the radial
  expression above after enforcing \(u(1)=1\).

Run with:

```sh
fpm run --example adaptive_helmholtz_bem_sphere
```
