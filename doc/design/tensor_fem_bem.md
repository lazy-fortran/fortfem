# Tensor FEM interior with a 3D BEM free boundary

`fortfem_tensor_fem_bem` couples the matrix-free tensor-product H1 operator of
the discretisation laboratory to the exterior Laplace problem through the
scalar Galerkin BEM of `fortfem_parametric_laplace_bem_3d`. It serves decision
criterion 1 (a general FEM interior plus a BEM free boundary). The interior is
any conforming curved hex mesh `tensor_mesh_t`, for example the 5-block
through-axis O-grid times phi torus, so the magnetic axis is an ordinary
interior line (criterion 2).

## Problem and coupled form

Let Omega be the meshed region with boundary Gamma and normal n out of Omega.
The interior field u solves `-div(K grad u) + c u = f`. The exterior field w
is harmonic outside Omega and decays at infinity. The transmission conditions
are `u - w = g_D` and `K grad u . n - grad w . n = g_N` on Gamma. With the
exterior energy map `S w = -d_n w`, the coupled weak form is

```text
a(u, v) + <S R u, R v> = (f, v) + <g_N, R v> + <S I g_D, R v>,
```

where R restricts the volume GLL face nodes to the BEM trace space and I
interpolates at the trace nodes. The discrete Steklov map is the existing
symmetric Schur complement

```text
S = W + (M/2 - K)^T V^{-1} (M/2 - K),
```

built by `build_exterior_dtn` without changes. The exterior energy map is
positive on nonzero constant traces, so constants are not in its kernel even
when the interior operator is a pure Neumann operator.

The Steklov Schur complement was chosen over the alternatives for the
present dense regime. The Johnson-Nedelec coupling needs only V and K, but it
is non-symmetric (GMRES) and keeps the flux unknown in the Krylov system. The
symmetric Costabel coupling keeps a symmetric indefinite saddle point with the
-V block, which needs MINRES and a two-sided preconditioner, but no dense
inverse. Eliminating the flux costs one dense LU of V per geometry,
O(N_flux^3), and gives an SPD system that runs PCG with the interior p-MG
unchanged. The same V, K, W and M blocks feed the symmetric Costabel form once
compressed operators are available, because an H-matrix or FMM does not supply
V^{-1} cheaply.

## Boundary trace map

`build_tensor_boundary_trace(mesh, trace, status)` is generic for 3D tensor
meshes. A cell face is on Gamma when its sorted four corner vertices occur in
exactly one cell. Each boundary face becomes one quadrilateral panel on the
reference square `[0,1]^2`, whose `(u,v)` tangent order makes `a x b` point out
of the cell for positive cell Jacobians. The panel geometry is the cell's
isoparametric geometry restricted to the face, so the BEM sees exactly the
discrete volume boundary. `panel_dof` maps face GLL nodes to trace unknowns, and
`trace_dof` maps trace unknowns to volume dofs, so R is a gather and R^T a
scatter-add. The builder rejects faces shared by more than two cells and
boundaries that are open or inconsistently oriented, and it returns the genus
from the Euler characteristic of the panel complex. The volume operator stays
matrix-free.

## BEM trace spaces on quadrilateral panels

`fortfem_parametric_bem_tensor` extends the existing parametric Laplace BEM
with tensor spaces instead of P1/P0 on triangles. It uses the same surface
contract `bem_surface_t%evaluate`, now on the square. It also keeps the same
generated frame, kernel and curl-density leaves, the same singular treatment
and the same `bem_operators_t` container:

- Potentials are continuous Q_p with Lagrange functions at the volume GLL
  nodes, so the BEM trace is the FEM trace and no interpolation is needed.
- Fluxes are discontinuous Q_pf at Gauss nodes, with `pf = p - 1` by default.
- Coincident pairs use a Duffy fan of four triangles centred at each target
  point, with radial coordinate `r = s^2`. Close pairs use the same fan about
  the closest source point. Far pairs use the tensor Gauss product rule. Close
  means a centre distance below `near_factor (d_i + d_j)`.
- The double layer integrates `d_nG (u_y - u_x)` and restores `-u_x/2`. The
  hypersingular form uses surface curls `b d_u phi - a d_v phi`.
- V and the K/M rows belong to the target panel's discontinuous flux unknowns,
  so OpenMP over target panels writes them without races. W uses atomic
  updates, as in the triangle path.

The representation formula (`bem_tensor_exterior_field`) evaluates the exterior
potential and gradient with the generated `bem_field_kernel`
(`tools/codegen/app/gen_tensor_fem_bem_products.f90`, which proves the
double-layer identity and the harmonicity of G). It uses regular quadrature
only and is valid a few panel sizes away from Gamma.

## Solver

`solve_tensor_fem_bem` runs PCG on A + R^T S R. The preconditioner is the
symmetric block factorisation for the interior I and trace B unknowns:

```text
y_I = MG(r_I);  z_B = (sigma S)^{-1} (r_B - A_BI y_I);  z_I = y_I - MG(A_IB z_B)
```

MG is one V-cycle of the existing Dirichlet p-multigrid with the trace nodes
as Dirichlet rows. The exact boundary Schur complement is `S + S_int`, where
`S_int` is the discrete interior Steklov-Poincare map. Both are H^(1/2) maps,
so it is replaced by `sigma S` with a dense Cholesky factor. Each iteration
costs two V-cycles, three matrix-free applies and O(N_trace^2) dense work.
`schur_scale` sigma defaults to 3, the minimum of the iteration counts scanned
on the O-grid torus (see `benchmark/lab/RESULTS.md`, FB1).

## API

```fortran
type(tensor_fem_bem_options_t) :: options      ! flux_degree, quadrature, schur_scale
type(tensor_fem_bem_t) :: fb
call build_tensor_fem_bem(mesh, coefficient, fb, status, options)
call tensor_fem_bem_rhs(fb, volume_load, rhs, status, neumann=g_N, dirichlet=g_D)
call solve_tensor_fem_bem(fb, rhs, tol, x, info)
call tensor_fem_bem_exterior_flux(fb, x(fb%trace%trace_dof) - gD_nodes, flux, status)
call bem_tensor_exterior_field(fb%trace%surface, fb%space, w_trace, flux, &
    points, potential, gradient, status)
call free_tensor_fem_bem(fb)
```

`volume_load` is a Neumann load, for example `tensor_load` on an operator
whose mesh has no Dirichlet rows. `g_N(x, n)` receives the panel normal, and
`g_D(x)` is evaluated at the trace nodes. `fb%steklov` is the dense S in trace
numbering and `fb%trace%trace_dof` gives R.

## Verification

- `test_lab_fem_bem_transmission`: exact transmission data from point charges
  outside the solid torus plus a non-harmonic quadratic (interior), and charges
  inside it (exterior). The data are exact on the discrete isoparametric
  domain, so the test has no geometry error. It checks volume L2 and H1 errors,
  the coupled Neumann error and a DtN-only Neumann error under p- and
  h-refinement, the trace topology (genus one, 4 n nphi panels) and the PCG
  iteration bound. Given arguments `n nphi p [pf nq nl sigma]`, it prints one
  scan row.
- `test_lab_fem_bem_current_sheet`: a vacuum-field free-boundary model. A
  prescribed divergence-free surface current `K = grad(Phi) x n` on Gamma
  produces `B = -grad(psi)` with psi the double-layer potential of Phi, i.e.
  `g_D = -Phi`, `g_N = 0`. The exterior B from the representation formula and
  the interior B on the axis core cells are compared with an independent
  Biot-Savart surface integral of K, whose own quadrature convergence is
  checked.

## Limits

- All BEM blocks are dense: O(N^2) storage and assembly, O(N_flux^3) LU of V
  and O(N_trace^3) Cholesky of the Schur block. Compression (ACA/H-matrices or
  FMM) is not implemented, and the measured crossover is recorded in RESULTS.
- The fan quadrature has no error bound. Convergence in the quadrature
  parameters must be demonstrated for each geometry, as the scans do.
- Only the single-valued scalar potential is coupled. A genus-one exterior also
  has a nonexact harmonic field (net toroidal or poloidal current), which needs
  the separate harmonic-sector data of `PARAMETRIC_LAPLACE_BEM_3D.md`.
- The preconditioner's `sigma S` block is spectrally equivalent to the exact
  Schur complement, but the constants depend on p and the mesh, so the
  iteration counts grow mildly with p.
