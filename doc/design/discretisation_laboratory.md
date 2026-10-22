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
and factored through FortFEM's native FortSparse hook with a symmetric hint selecting native LDLT. A
procedure argument can replace the coarse solve, for example with AMG. PCG
confirms convergence using a fresh operator residual. Nonzero Dirichlet values
are lifted before solving the masked interior system.

B-spline and Fourier builders plug into the line contract. The present
isoparametric mesh and p-coarsening adapter use GLL lines, so mixed-family
operator construction still needs an adapter selecting per-direction tables
and transfer rules. The shared contract itself requires no changes.

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
