# Discretisation laboratory

## Decision criteria

1. A clean, general FEM interior plus BEM free boundary, using FortFEM's 3D
   Laplace BEM / Dirichlet-to-Neumann coupling.
2. No magnetic-axis problem: axis-regular bases and charts, such as Cartesian
   or cylindrical meshes through the axis (multi-block O-grid or unstructured
   curved quads times phi, or Cartesian hexes). Polar singular treatment is
   not the structural approach. Polar `(s,theta)` charts are comparison cases.
3. Field alignment is optional, an accelerator rather than a structural
   requirement.
4. Time to solution at matched accuracy must reach or beat the fastest
   VMEC-class codes (VMEC++, GVEC) on nested cases while also handling
   non-nested and island cases. KIN6D must not clone VMEC.

## Free boundary

The FEM interior couples to the exterior Laplace problem through the 3D BEM
Dirichlet-to-Neumann Schur complement on the isoparametric boundary of any
tensor hex mesh. Its design, API and limits are in
[tensor_fem_bem.md](tensor_fem_bem.md).

## Abstraction

`tensor_line_tables_t` is the common, unchanged line-table contract. It carries
basis values and parameter derivatives at quadrature points, weights including
cell length, and periodic or nonperiodic dof identification. GLL Lagrange,
B-spline and Fourier builders supply this contract. A Fourier line represents
one global cell and is suitable only along an exact symmetry direction.

`tensor_mesh_t` separates field degree from geometry degree. The geometry is
an array of physical coordinates at reference GLL nodes. Cell corners use
binary tensor ordering, with the first direction varying fastest. An explicit
cell-to-global-dof map supports arbitrary conforming quads and hexes, including
periodic identifications. Without an explicit map, sorted global corner IDs
and integer tensor barycentric weights identify shared vertex, edge and face
nodes under reference rotations and reversals. This constructor assumes a
regular conforming cell complex, with a face identified by its corner set.
Periodic multigraphs with parallel edges need the explicit gather map, as used
by the structured patch generator. Boundary masks are caller-owned for general
meshes. Patch and O-grid generators supply their boundary masks.

The operator gathers, contracts one direction at a time, applies a generated
pointwise metric and reaction kernel, contracts test functions, and scatters
by cell colours. The complexity per cell is O(p^(d+1)). The exact diagonal uses
squared and mixed line tables with the same tensor contractions. A separately
scheduled, unfactored quadrature assembly provides the CSR oracle.

Precomputed metrics trade storage for arithmetic. The recompute mode retains
Jacobians and sampled physical coefficients and regenerates the metric during
apply, without calling user procedures in the parallel hot loop. Both modes
use the same isoparametric geometry. Recomputing the full geometry interpolation
is a separate optimisation and is not what this mode measures.

The p-multigrid path rediscretises at `p -> ceil(p/2) -> ... -> 1`, retaining
geometry degree. GLL interpolation is tensorised. Each fine global dof has
one owner cell, making restriction exactly the transpose of global
prolongation. Jacobi or Chebyshev(Jacobi) smoothing uses a short Lanczos
spectral estimate with a safety factor. The p=1 CSR matrix is converted to CSC
and factored through FortFEM's native FortSparse hook with the SPD hint, which
selects native Cholesky. The masked coarse matrix is SPD, and the pivoted
indefinite LDLT used before took 42 s on a 5,185-dof 2D coarse level. A
procedure argument can replace the coarse solve, for example with AMG. PCG
confirms convergence using a fresh operator residual. Nonzero Dirichlet values
are lifted before solving the masked interior system.

B-spline and Fourier builders plug into the same line contract and the same
operator; see "Mixed lines in the MF operator" below.

## MF core status

The CPU source includes a GLL builder, structured and general conforming maps,
a five-block axis-crossing disk and shaped periodic torus generator, curved
geometry, precomputed and recomputed diffusion metrics, coloured OpenMP apply,
exact diagonal, CSR oracle, PCG and p-multigrid. The torus uses nfp=5 smooth
Cartesian deformation through its axis. Geometry callbacks build structured
patch coordinates, and general cells accept coordinate arrays.

`solve_tensor_diffusion(mesh, coefficient, rhs, tol, x, info, bc=...)` accepts
an assembled load vector and optional nodal boundary values. The coefficient
callback returns symmetric physical K and reaction c. `tensor_load` forms the
load vector from a scalar forcing callback. Fully periodic Poisson needs a
caller-selected gauge or positive reaction before PCG can solve an SPD system.

Named native tests are `test_lab_gll_line`, `test_lab_tensor_mesh`,
`test_lab_tensor_apply`, `test_lab_tensor_solver` and `test_lab_tensor_mms`.
Measured evidence and remaining acceptance gates belong in
`benchmark/lab/RESULTS.md`; implementation alone does not establish convergence,
scaling or a VMEC-class time-to-solution claim. Interior/BEM coupling is a
laboratory decision criterion, not an implemented boundary operator here.

The validated periodic meshes have no dof self-identifications within a cell.
One-cell periodic self-identifications are outside the exact-diagonal test
scope. Coarsening of explicit periodic multigraph gather maps also needs
separate validation because the automatic coarse constructor uses corner
entities of a regular cell complex.

## Spline and Fourier lines

`fortfem_bspline_lines` and `fortfem_fourier_lines` fill the shared
`tensor_line_tables_t` contract. Quadrature weights include span lengths and
`d` differentiates in the parameter coordinate. These are polynomial spline
and Fourier tables, with geometry and rational NURBS weights supplied by the
operator layer. Construction is setup work, outside the operator hot path.

`build_bspline_line(degree, line, status, ...)` accepts exactly one of
`nelem` (uniform spans on [0,1]) and `knots`. Open `knots` is a complete,
nondecreasing, clamped knot vector, including endpoint multiplicities p+1.
Interior multiplicities determine continuity and may range from 1 to p.
Periodic `knots` instead gives strictly increasing span breaks, including
both endpoints. `continuity` defaults to p-1 for generated knots and periodic
breaks and accepts 0 through p-1. If supplied for an open knot vector it
checks its interior multiplicities. Degree zero builds piecewise constants.
`nq` defaults to p+1 Gauss points per positive span. Invalid inputs return a
nonzero integer status. Nonuniform spans are supported in both cases.

The builder calls the existing `evaluate_bspline_basis` Cox-de Boor routine.
Periodic knots extend by translation, and local support maps to global dofs
modulo the number of periodic functions. Several local functions can map to
the same dof on a small periodic mesh, so consumers must scatter-add them.

`build_bspline_line_pair` requires p >= 1 and returns V0, V1 and a dense
coefficient incidence matrix D with shape (V1%ndof,V0%ndof). V1 removes the
first and last knot, lowers degree to p-1 and continuity by one (C0 becomes
discontinuous), and scales its functions as
`M_j = p B_{j,p-1}/(t_{j+p+1}-t_{j+1})`, with indices in the V0 knot vector.
Thus `(D u)_j = u_{j+1}-u_j`, wrapping at the periodic seam, and the
reconstructed derivative equals the V1 field at every point inside a span.
Open pairs have one fewer V1 dof. Periodic pairs retain equal dimensions and
a constant cokernel, as required by the topology of a circle. Tensor products
of the directional V0/V1 pairs provide the patch H1/H(curl)/H(div)/L2 spaces.

`build_fourier_line(M, line, status, nq, nfp)` uses one periodic cell on
[0,2*pi/nfp), ordered `1, cos(nfp*x), sin(nfp*x), ...`. Both nfp and nq are
optional and default to 1 and 2M+1, respectively. The equally weighted
trapezoidal rule requires nq >= 2M+1 and exactly integrates retained basis
products. Evaluation and derivatives reuse FortFEM's Fourier mode registry.
The mass diagonal is the period for the constant and half the period for
other modes, and stiffness divided by mass has eigenvalues (k*nfp)^2.
`build_fourier_line_pair` chooses V1 to be the derivative image, omitting
the constant. Each cos/sin block of D maps to `(-k*nfp*sin, k*nfp*cos)`.
For M=0, V1 has zero functions. This zero-mean V1 convention omits the harmonic
one-form of a full periodic de Rham complex. Consumers needing that harmonic
component should retain the constant using `build_fourier_line` for V1.

GLL derivative/histopolation pairs are not supplied here. The MF1 nodal line
builder provides V0, but a compatible degree-reduced edge basis remains a
separate extension.

`test_lab_spline_lines` checks partitions, polynomial moments, all basis
commuting identities, periodic wrapping and cardinal values, Fourier products
and eigenvalues, and 1D and 2D mass/stiffness against the existing sparse
B-spline assembly on an affine nonuniform patch, with a 1e-13 assembly tolerance.

## Mixed lines in the MF operator

`fortfem_line_mesh` builds meshes from one line table per direction and a
per-cell 1D cell index, so GLL, B-spline (any degree and continuity, open or
periodic) and Fourier lines and their products use one gather, one
sum-factorised interpolation, the same generated metric kernel and one
scatter. Cells are knot spans for splines and one global cell for Fourier.
Geometry is stored per cell as Lagrange samples of degree q in every
multi-cell direction. In a single-cell Fourier direction it is stored at the
line's nq equispaced quadrature points, with the periodic spectral
differentiation matrix, so the geometry does not depend on the number of
retained modes and is identical on mode-coarsened levels.
`line_ogrid_mesh` glues five open-knot patches through the axis with C0
interfaces through the general gather map, using MF1's O-grid layout.
`tensor_boundary_projection` supplies Dirichlet data for non-nodal bases by a
face-wise L2 projection.

`fortfem_line_hierarchy` builds nested hierarchies. B-spline levels halve the
patch spans by knot insertion, and periodic spline phi lines halve their spans
while at least degree+1 spans remain. Fourier levels truncate modes, and GLL
levels lower the degree p -> (p+1)/2 on the same cells. A single B-spline box
supports per-direction semi-coarsening for aligned anisotropy.
`line_mesh_prolongation` forms the global prolongation from the 1D matrices,
and `build_line_multigrid` reuses MF1's Chebyshev-Jacobi smoothing, its
transposed restriction and its coarse LDLT. The Galerkin identity
P^T A_fine P = A_coarse holds to roundoff for O-grid spline and Fourier
hierarchies when geometry and quadrature agree between the levels.

C1 polar splines (`fortfem_bspline_polar`, Toshniwal-Hughes extraction) enter
only as a comparison for criterion 2. The driver assembles the tensor
operator on the polar map, restricts it by the extraction and solves it
directly. `fortfem_lab_igaax` writes one CSV row per run, and
`test_lab_line_operator_slow` and `test_lab_line_hierarchy` are the named tests.


The IGA-AX section of `benchmark/lab/RESULTS.md` holds the evidence. At
matched L1 error, GLL with p-multigrid was the fastest family (1.26 s for
rel L2 1.8e-9 at 16 threads). The O-grid through the axis kept rate p+1
without any near-axis degradation for both GLL and splines. Spline
h-multigrid with Jacobi-type smoothing is h-robust but not p-robust.

## Nonmatching periodic seams

A mesh whose toroidal lines follow field lines closes a period with a poloidal
twist that is in general not a multiple of the cell size, so its two periodic
faces do not match node for node. Alignment stays an optional accelerator
(criterion 3), and the seam closure exists so that it can be measured.
`close_tensor_seam` in `fortfem_tensor_seam` takes the two face lists, a rigid
seam motion and optional lattice images. A matching seam is closed by DOF
identification, bit for bit the ordinary periodic gather map. Otherwise the
slave face DOFs are eliminated by an affine constraint `u_s = C u_m`, and every
operator acts as `P^T A P` with `P = [I; C]`, matrix free through
`apply_seam_operator`. The assembled `assemble_seam_operator` is the oracle and
the coarse matrix.

Two constraints are available. `SEAM_INTERPOLATION` sets slave nodes to the
master trace at their images, as flux-coordinate-independent codes do along
field lines ([Hariri and Ottaviani 2013](https://doi.org/10.1016/j.cpc.2013.06.005);
[Stegmeir et al. 2016](https://doi.org/10.1016/j.cpc.2015.09.016)). The master
trace has kinks at master cell edges inside slave cells, so the closed space is
not continuous across the seam and the pointwise condition is not optimal for
p >= 2. The default `SEAM_MORTAR` is the mortar condition of Bernardi, Maday and
Patera (1994): slave test functions, lowered by one degree along a Dirichlet
rim, are orthogonal to `u_slave - u_master o T`
([Ben Belgacem 1999](https://doi.org/10.1007/s002110050468) analyses the
multiplier form). It is eliminated as a constraint, so it costs the same per
apply as interpolation and needs no saddle-point solver. The flux-tube
twist-shift condition ([Beer, Cowley and Hammett 1995](https://doi.org/10.1063/1.871232))
states the same identification for a Fourier binormal coordinate, whereas a C0
GLL space needs it in real space.

The overlap integrals are exact up to quadrature on curved cells. Master faces
are split into nearly bilinear parts, whose charts invert reliably near their
squares, and slave faces into still flatter parts. A slave part that is larger
than an eighth of a master part it meets is halved against that master part
only; at a quarter, fan nodes of a rigid twist fell outside the slave part.
Each slave part is clipped in the chart of a master part, and a Gauss rule on
the fan of curved triangles from the polygon centre integrates the overlap, so
no master or slave kink lies inside a quadrature cell. Clipping a large, nearly
affine slave face in a curved master chart, as a shear-free twist produces, met
spurious fold preimages, and the size bound removes that case. A global bound,
a quarter of the smallest master part for every slave part, built 45 million
overlap nodes for the five slave faces of the p=2 L3 seam and took 183 s; the
pairwise bound keeps the refinement where the small master parts are. Each fan
triangle uses the p+q+2 Gauss points per direction of the slave mass rule,
where q is the geometry degree; the earlier 2p+q+4 points changed the
reproduction and operator checks only at roundoff. The fan keeps its nodes
inside the nearly convex overlap, where both charts invert reliably; a node
outside either cell fails the setup instead of entering the constraint. The
discrete rims of the two faces differ by the geometry error of the curved
boundary. A master cell on a Dirichlet rim is therefore not clipped on its rim
side, and slave slivers beyond it see the polynomial extension of its trace.
Without this, constants are not reproduced next to a curved rim. The overlap
cubature depends only on the geometry, which all p-levels share, so
`coarsen_tensor_seam` reuses it and p-multigrid restricts and prolongs on
consistently constrained coarse spaces. The setup is dense in the slave rows of
the seam, O(n_slave^3), which is acceptable for one cross-section but not for a
seam with many thousands of slave DOFs.

L3 evidence ([V], sCluster job 3356493, `benchmark/lab/RESULTS.md` on the L3
integration branch): on the exact unit disk with n=1 and 20 phi cells, the
aligned O-grid with this mortar seam is not more accurate than the periodic
O-grid at matched p, except at p=2 and ratio 1e3, and is seven times less
accurate at p=6 and ratio 1e3 (3.5e-7 against 4.8e-8). It needs more PCG
iterations, and its wall time is 106 to 362 times larger because seam setup
dominates (9, 36 and 111 s at p=2, 4 and 6 on four threads, before the
robustness fix above). Already at ratio 1 the aligned error is 1.3 to 4.4
times the disk error. Field alignment through a nonmatching seam is therefore
not an accelerator on this case, which is consistent with criterion 3 keeping
alignment optional. The seam stays available for charts that need it, such as
KIN6D field-line charts with irrational iota, and its setup cost is the open
item (overlap cubature on finely split parts, and the dense slave solve).
