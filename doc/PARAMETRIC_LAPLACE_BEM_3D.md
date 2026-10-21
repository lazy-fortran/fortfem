# General parametric-surface Laplace BEM

[V] `fortfem_parametric_laplace_bem_3d` provides a scalar Galerkin path for a
fixed, connected, closed, oriented triangular surface. It does not assume
axisymmetry. Potential traces are continuous P1 and normal fluxes are P0.
Geometry is evaluated through a caller-owned map and may be curved to any
order, including an exact analytic map. High-order **trace spaces** are not
implemented by this API.

## Surface and operator API

[V] Extend `bem_surface_t` and implement the pure procedure
`evaluate(self,panel,uv,x,a,b,status)`. The reference triangle is
`uv(1)>=0, uv(2)>=0, sum(uv)<=1`; `a` and `b` are the position derivatives
with respect to its two coordinates. The map must be regular, agree between
neighboring panels, and orient `a cross b` out of the bounded region. Pure
maps must be safe to call concurrently. `faces(3,npanels)` supplies the
matching P1 vertex numbering and oriented topological connectivity.

```fortran
call assemble_parametric_laplace_bem(surface, faces, degree, op, status)
call build_exterior_dtn(op, exterior, status)
call exterior%apply(dirichlet_trace, outward_flux, exterior_load, status)
call exterior%matrix(steklov, status)
```

[V] Assemble and build once per fixed geometry, connectivity and quadrature
choice. Each later apply is a matrix action, without assembly or factorization.
The container has `v`, `k`, `adjoint`, `w`, `mass`, `area` and `genus`.
`v` acts on panel flux coefficients; `k` and `mass` map vertex traces to
panel tests; `w` maps vertex traces to vertex tests. `adjoint=transpose(k)`.
All operators and accepted DtN actions use fp64.

[V] With `n` pointing into the exterior and `G=1/(4*pi*|x-y|)`, the physical
outward derivative satisfies `V q=(K-M/2)u`. The symmetric exterior energy map
is `S=W+(M/2-K)^T V^{-1}(M/2-K)`. `apply` returns this positive energy load
`S u` as well as the separate physical flux `q`. These have different discrete
spaces and must not be interchanged. `matrix` returns a copy suitable for the
FEM restriction `R^T S R`. A caller must add its known exterior drive using
the same normal and trace conventions.

[V] Status 1 rejects inputs or topology; status 2 rejects geometry or
nonfinite assembly values. The DtN builder additionally propagates the dense
solver status. An unsuccessful output must not be applied. Changing a public
operator array requires rebuilding the cached DtN object.

## Quadrature and verification limits

[V] Coincident-panel integration uses a source Duffy fan centered at each
target point. Close-panel integration uses a source fan about a projected
reference point and a squared radial coordinate to cluster samples. The
double layer integrates `d_n G (u(y)-u(x))` and restores the exact
closed-surface term `-u(x)/2`. The hypersingular form uses surface curls,
reducing it to a weakly singular integral. Far pairs use product triangle quadrature. Assembly is OpenMP;
shared hypersingular entries use atomic additions.

[V] `degree` controls quadrature, not approximation order. This implementation
has no tolerance-controlled near-pair error bound or Sauter-Schwab edge/vertex
pair rule. A user must demonstrate quadrature convergence and mesh convergence
for the actual shell, particularly for close but nonadjacent panels.

[V] FortSym generates the Green, normal, Jacobian, curl-density and integrand
leaves, including JVP leaves, from `tools/codegen/app/gen_parametric_bem_products.f90`.
Generated output is checked in. The assembled operator currently has no
geometry-JVP/VJP API; fixed near/far selection and the surface-map derivatives
must be included before claiming differentiable assembly.

[V] `parametric_bem_residual(op,u,q,panel_indicator,status)` reports the
area-scaled discrete first-Calderon equation defect. With the assembly operator
used to solve for `q`, it measures algebraic error only. With an independently
more accurate quadrature operator on the same mesh, it can reveal quadrature
error. It is neither a continuous error bound nor a discretization certificate;
a refined-mesh estimator remains necessary. The independent sphere and torus
oracles are in `test/test_parametric*`.

[V] Manufactured flux tests compare the P0 coefficients with the exact
area-averaged normal derivative on each panel. Unresolved variation inside a
panel is excluded from that coefficient error. A continuous field-error norm
and estimator are still needed for a total discretization-error claim.

## Topology and remaining integration

[V] Connectivity checks reject open edges, nonmanifold edges, inconsistent
edge orientation and disconnected surfaces, and report the orientable genus.
The scalar hypersingular form annihilates the constant trace; the decaying
exterior Dirichlet map itself does not require removing that constant.

[V] A single-valued scalar potential has zero periods. A genus-one exterior
can also carry a nonexact harmonic one-form. Its circulation, and the matching
interior harmonic flux data, are separate data, not extra null vectors of the
scalar single-layer matrix. This API does not construct or couple those fields.
FortFEM's existing `cell_complex` and `harmonic_period_normalization` modules
provide topology and period-normalization primitives for a subsequent adapter.

[V] This API is not a complete general magnetic exterior closure. High-order
trace spaces, a reliable refined-mesh estimator, harmonic-sector coupling,
assembled geometry derivatives, compressed assembly/solves and controlled mixed
precision remain open. The existing per-mode axisymmetric API remains the
consumer's established route until these gaps and its accuracy gates close.

[V] Operator and weak-form conventions follow the
[Bempp boundary-operator definitions](https://bempp.com/handbook/api/boundary_operators.html)
and the existing FortFEM `toroidal_mode_boundary_operators_2d` sign conventions.
