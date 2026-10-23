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
only; at a quarter, fan nodes of a rigid twist fell outside the slave part. Each slave part is clipped in the chart of a master part, and a Gauss
rule on the fan of curved triangles from the polygon centre integrates the
overlap, so no master or slave kink lies inside a quadrature cell. Clipping a
large, nearly affine slave face in a curved master chart, as a shear-free
twist produces, met spurious fold preimages, and the size bound removes that
case. A global bound, a quarter of the smallest master part for every slave
part, built 45 million overlap nodes for the five slave faces of the p=2 L3
seam and took 183 s; the pairwise bound keeps the refinement where the small
master parts are. Each fan triangle uses the p+q+2 Gauss points per direction
of the slave mass rule, where q is the geometry degree; the earlier 2p+q+4
points changed the reproduction and operator checks only at roundoff. The fan keeps its nodes inside the nearly convex overlap, where both charts invert reliably; a node outside either
cell fails the setup instead of entering the constraint. The discrete rims of the two faces differ by
the geometry error of the curved boundary. A master cell on a Dirichlet rim is
therefore not clipped on its rim side, and slave slivers beyond it see the
polynomial extension of its trace. Without this, constants are not reproduced
next to a curved rim. The overlap cubature depends only on the geometry, which all p-levels
share, so `coarsen_tensor_seam` reuses it and p-multigrid restricts and
prolongs on consistently constrained coarse spaces. The setup is dense in the
slave rows of the seam, O(n_slave^3), which is acceptable for one cross-section
but not for a seam with many thousands of slave DOFs.
