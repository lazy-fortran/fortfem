---
title: Differentiation strategy
---

# Differentiation strategy

FortFEM uses the same public derivative contract as FortNum:

```fortran
call foo(...)
call foo_jvp(..., x_dot, ..., y_dot)
call foo_vjp(..., y_bar, ..., x_bar)
```

The unsuffixed routine evaluates the primal operation. A JVP propagates one
tangent direction without forming a Jacobian. A VJP propagates one output
cotangent to all active inputs without forming a transpose Jacobian. Output
bars are returned values rather than hidden global tape state.

## Method selection

Derivative implementations are selected per kernel, not per application:

1. FortSym generates compact analytical JVPs and VJPs for fixed algebraic
   geometry, basis, quadrature, and local-form kernels. These compile with GCC
   and are the default production path.
2. Matrix and state solves use implicit differentiation. For
   \(A x=b\),
   \(A\dot x=\dot b-\dot A x\), while one transpose solve gives
   \(\bar b=A^{-T}\bar x\) and
   \(\bar A_{ij}=-\bar b_i x_j\) on the active sparse pattern. FortFEM does not
   differentiate through a third-party factorization.
3. Enzyme is an optional implementation oracle and performance candidate for
   kernels that contain control flow or become impractical to generate.
   Enzyme is not required by GCC users.
4. Central differences, dot identities, manufactured solutions, and external
   discretization oracles test behavior independently. Agreement between
   Enzyme and FortSym alone is not accepted as a test oracle.

Forward mode is preferred when the number of input directions is small.
Reverse mode is preferred for scalar or low-dimensional objectives with many
parameters, including NURBS control points, rational weights, and mesh-motion
degrees of freedom. Hybrid paths compose analytical local products, sparse
implicit adjoints, and application-level reverse accumulation.

Fixed-mask Dirichlet elimination has matching analytical products. The
constraint mask is a discrete, inactive topology choice; CSC values, the
right-hand side, and prescribed values are active. Both products factor only
the reduced free system. Reverse mode scatters free--free and
free--constrained matrix cotangents back into the original CSC pattern and
returns the direct plus elimination-mediated boundary-value cotangent.

## H1 geometry tournament

The optional CMake test
`fortfem_enzyme_h1_geometry_products` differentiates the same local
isogeometric H1 contribution with FortSym and Enzyme. It checks:

- the analytical JVP against a two-sided directional difference;
- the analytical JVP/VJP dot identity;
- Enzyme forward and reverse products against the analytical products; and
- repeated-call timings on an identical, deterministically varying geometry.

Configure it with:

```sh
cmake -S . -B build-enzyme -G Ninja \
  -DFORTFEM_ENABLE_ENZYME=ON \
  -DFORTFEM_ENZYME_REQUIRED=ON \
  -DFORTNUM_ENZYME_PLUGIN=/path/to/LLVMEnzyme-22.so
cmake --build build-enzyme \
  --target fortfem_enzyme_h1_geometry_products_build
ctest --test-dir build-enzyme --output-on-failure -R fortfem_enzyme
```

The fixture reuses FortNum's revision-pinned Flang/LLVM/Enzyme pipeline.
Continuous integration builds a pinned Enzyme plugin and requires the real
forward and reverse test to pass.

On the reference LLVM 22 toolchain, a one-million-call local tournament found
Enzyme forward mode comparable to the generated JVP and Enzyme reverse mode
about 1.3 times the generated VJP. These are toolchain-specific observations,
not permanent dispatch thresholds; CI reports fresh values on every run.

## Current end-to-end path

The `iga_shape_sensitivity` example composes NURBS geometry products, sparse H1
assembly products, and exact implicit sparse-solve products. One reverse solve
returns the gradient with respect to every control point and weight. A complete
reassemble-and-refactor finite difference independently checks the resulting
geometry-to-state adjoint.

Cartesian PML coefficient maps expose analytical JVPs and VJPs for both scalar
Helmholtz and Maxwell curl--curl forms. FortSym generates the complex local
products. Their VJPs use the real-vector complex convention
\(\operatorname{Re}(\bar y^H\dot y)
=\operatorname{Re}(\bar x^H\dot x)\), so reverse mode applies the conjugate
transpose required by optimization over real and imaginary design components.
The preceding element-stretch map is differentiable with respect to every mesh
vertex, physical and outer boundary coordinate, wave number, and attenuation
amplitude. Its reverse product accumulates contributions at shared vertices.
At the physical/PML interface the active-set branch is intentionally
piecewise: interior cells have zero stretch derivative, while derivatives in a
PML layer are valid as long as a perturbation does not cross the interface.

The planar Helmholtz DtN trace operator has analytical trace, wave-number, and
period JVPs and VJPs. FortSym generates the propagating and evanescent modal
square-root products, while FortNum's FFT adjoint supplies the dense trace
reverse product without forming the Fourier matrix. A mode exactly at cutoff
is rejected by parameter products because the square-root derivative is
singular there; the primal operator remains defined.
The corresponding fixed-mesh P1 FEM--DtN state solve composes those products
with the complex implicit FortSparse adjoint. It differentiates wave number,
period, complex volume loads, and complex Dirichlet data. Reverse elimination
accounts for the dependence of the reduced right-hand side on both the
pre-elimination matrix and prescribed boundary values.

The biperiodic planar Maxwell capacity operator follows the same public
`foo`, `foo_jvp`, and `foo_vjp` contract. Its products differentiate both
tangential trace components, the wave number, and the two cell lengths. The
FFT remains matrix-free, the modal square-root derivative is shared with the
FortSym-generated Helmholtz kernels, and the two-by-two TE/TM coupling is
propagated analytically. The assembled weak boundary form composes these
products with the cell-area quadrature weight, including both length
derivatives. Reverse products use the real inner product
`real(sum(conjg(output_bar)*output_dot))`; exact cutoff modes are rejected
because their parameter derivative is singular.

The Nedelec trace coupling exposes the smooth pullback
`sampling^T * capacity * sampling` as a separate differentiable primitive.
Its JVP accepts a sampling perturbation, and its VJP returns a dense
sampling-matrix cotangent together with the Maxwell parameters. This keeps
the reciprocal bilinear transpose used by the weak form distinct from the
conjugate transpose required by the real-complex reverse inner product.
Mesh point-location and boundary-face selection remain explicit fixed-topology
operations; geometry products can be composed on either side of the sampling
cotangent without pretending those discrete decisions are differentiable.

Spherical scalar Helmholtz and Maxwell DtN modes expose eigenvalue and
applied-operator JVPs and VJPs for modal traces, wave number, and sphere
radius. The derivative of the outgoing Hankel logarithmic derivative is
evaluated analytically from the spherical Bessel differential equation and
shared by the scalar and TE/TM maps. This avoids differentiating the
special-function recurrence, keeps the implementation available with GCC,
and preserves the exact TE/TM impedance duality in the primal map.

A discrete planar Maxwell FEM--DtN state primitive composes an assembled
curl--curl volume matrix, the differentiable Nedelec trace pullback, and the
complex FortSparse implicit solve. Its JVP accepts volume-matrix, trace-map,
load, wave-number, and cell-size perturbations. One adjoint solve returns
cotangents for every volume-matrix entry, every trace-sampling entry, the
complex load, and all three scalar boundary parameters. This algebraic
boundary lets element and moving-mesh geometry products compose with the
state adjoint without making mesh topology or point location implicit.

The tetrahedral Nedelec covariant Piola map provides products for the
Jacobian, reference values, and reference curls. Its implementation composes
FortNum's FortSym-generated determinant and inverse JVP/VJP kernels with the
closed covariant value and curl transformations. Consequently element and
mesh geometry adjoints remain analytical and GCC-compatible while sharing
the same guarded linear-algebra primitives used elsewhere.

The tetrahedral RT contravariant Piola map likewise exposes Jacobian,
reference-value, and reference-divergence products. It composes the
FortSym-generated determinant products with the closed H(div) value and
divergence transformations, providing the geometry primitive required by
mixed Poisson and flux-conservative Ampere discretizations.

Arbitrary-order tetrahedral RT div--div-plus-mass element matrices compose
that map with quadrature. Their JVPs cover all twelve vertex coordinates and
both operator coefficients. Their VJPs use one reverse quadrature sweep to
return the complete geometry and coefficient cotangents, making reverse mode
independent of the number of active mesh coordinates.
The global CSC products preserve the merged sparse pattern, reverse RT face
orientation and permutation transforms, and accumulate shared-vertex
cotangents across adjacent tetrahedra.

The algebraic RT--DG mixed state solve exposes flux-mass, divergence, and load
JVPs and VJPs. It assembles the saddle system
`[M, -B^T; B, 0]`, reuses one primal factorization for the tangent solve, and
uses one transposed factorization for the complete reverse product. The VJP
returns cotangents aligned with the input CSC patterns, so assembly products
compose without converting sparse design variables to dense matrices.
For affine RT elements the divergence block has no mesh derivative:
`det(J) div(J v_hat/det(J)) = div(v_hat)` exactly. Moving-mesh mixed Poisson
adjoints therefore propagate geometry through the RT mass block and physical
load, while the algebraic API still permits an independently active
divergence operator for non-affine or externally assembled discretizations.
The tetrahedral mixed-Poisson state wrapper performs that composition
end-to-end. Its JVP accepts every mesh-coordinate direction, the flux-mass
coefficient direction, and an assembled DG load direction. Its VJP returns
all shared vertex cotangents, the material cotangent, and the load cotangent
using a single sparse adjoint solve followed by one RT assembly reverse sweep.
Keeping the load explicit avoids silently assuming derivatives for an opaque
Fortran source callback; spatially differentiable source representations can
compose at this boundary.
Sampled tetrahedral RT vector loads provide that representation. Their
products include Eulerian vector-source gradients, contravariant Piola
geometry, determinant quadrature, arbitrary-order face transforms, and
shared mesh accumulation. The VJP returns source-sample and mesh-coordinate
cotangents without differentiating an opaque callback.

Arbitrary-order tetrahedral Nedelec curl--curl-plus-mass element matrices
compose those Piola products with quadrature. Their JVPs cover all twelve
vertex coordinates and both material coefficients. A single analytical
reverse quadrature sweep returns the corresponding vertex and coefficient
cotangents, rather than requiring one forward sweep per design coordinate.
This is the preferred element path for large moving-mesh design spaces.
Sampled Nedelec vector loads use the same Eulerian source contract as RT
loads. Their products reverse covariant Piola values, determinant quadrature,
arbitrary-order edge and face transforms, and source point motion. This
closes the differentiable forcing path needed by curl--curl and Ampere state
adjoints.
The constrained sampled-state wrapper composes that forcing with global
curl--curl-plus-mass assembly and fixed-mask elimination. One sparse adjoint
solve returns mesh, material, Eulerian source, and prescribed tangential-DOF
cotangents for Ampere and Maxwell optimization.
Mixed Poisson has the corresponding sampled discontinuous scalar-load
contract. Its analytical products include determinant quadrature, physical
quadrature-point motion, arbitrary-order DG moments, and source parameters.
The RT--DG sampled-state wrapper composes those products with the saddle
operator and one transpose solve, returning mesh, flux-mass coefficient, and
physical-source cotangents without differentiating a sparse factorization.
The two-dimensional FEEC foundation uses the same contract: covariant
Nedelec and contravariant RT triangle Piola maps expose analytical JVPs and
VJPs. Their determinant and inverse products reuse FortNum's
FortSym-generated kernels, so GCC builds retain the fast analytical path
without requiring an autodiff compiler.
Arbitrary-order triangle RT div--div-plus-mass elements compose those map
products with determinant quadrature and material coefficients. Their reverse
product accumulates both vertex coordinates and coefficients in one element
sweep, preparing the global two-dimensional H(div) assembly for shape
adjoints.
The global triangle RT CSC products preserve the orientation-aware merged
sparsity pattern. Reverse assembly looks up only participating CSC entries
rather than allocating a dense global adjoint, then accumulates shared-vertex
and material cotangents through the element VJP.
Arbitrary-order triangle Nedelec curl--curl-plus-mass elements likewise
differentiate vertex geometry, the curl coefficient, and every entry of an
anisotropic mass tensor. The tensor product is intentionally explicit so
material and curvilinear-coordinate optimization do not collapse to an
isotropic scalar parameter.
The global triangle Nedelec CSC products preserve edge orientations and the
merged sparse pattern while accumulating shared mesh vertices, the curl
coefficient, and the full anisotropic mass tensor. As for RT, reverse
assembly queries participating CSC entries and does not materialize a dense
global cotangent.
Triangle RT sampled vector loads follow the Eulerian source convention used
in three dimensions. Their products include physical quadrature-point motion,
contravariant Piola values, determinant weights, edge orientation, source
parameters, and shared mesh vertices.
Triangle Nedelec sampled loads expose the identical Eulerian source API with
covariant Piola products. This makes physical vector forcing interchangeable
between the two-dimensional H(div) and H(curl) state wrappers while retaining
the correct geometry pullback for each space.
The constrained triangle Nedelec sampled-state wrapper composes this forcing
with anisotropic curl--curl assembly and the differentiable sparse solve.
One transpose solve returns cotangents for shared mesh vertices, the curl
coefficient, the full mass tensor, physical source samples, and prescribed
tangential degrees of freedom.
The constrained triangle RT sampled-state wrapper provides the matching
H(div) path. It differentiates the SPD div--div-plus-mass system with respect
to shared mesh vertices, both material coefficients, physical source samples,
and prescribed normal-flux degrees of freedom.
The full polynomial triangle sequence is covered at element level as well:
second-kind Nedelec and BDM curl/div-plus-mass elements share one analytical
product implementation, selecting the appropriate covariant or
contravariant Piola pullback. This avoids duplicating the quadrature reverse
sweep while exposing family-specific public JVP/VJP names.
Their global CSC products also share one orientation-aware implementation.
The reverse sweep queries only participating sparse entries and accumulates
shared vertices and both coefficients for BDM and second-kind Nedelec without
a dense global cotangent.
BDM and second-kind Nedelec sampled vector loads also share one Eulerian
forcing implementation. The family switch changes only the Piola map; source
motion, determinant quadrature, orientation transforms, source cotangents,
and shared-vertex accumulation are common.
Their constrained sampled-state wrappers likewise share one end-to-end
primal/JVP/VJP pipeline. The products compose full-vector CSC assembly,
Eulerian sampled forcing, and the implicit sparse-solve adjoint, returning
cotangents for every shared mesh vertex, both scalar coefficients, source
samples, and prescribed boundary degrees of freedom without forming a state
Jacobian.
Fixed-reference-point RT, BDM, and first- and second-kind Nedelec observations
complete the state-to-objective path. Their shared products differentiate the
physical vector value and divergence or curl with respect to element vertices
and every local coefficient. They reuse the corresponding Piola products;
inversion of a fixed physical observation point is deliberately kept as a
separate geometry operation.
The affine triangle physical-to-reference inverse now has its own products.
They differentiate all six vertex coordinates and both physical point
coordinates by composing FortNum's FortSym-generated two-by-two inverse
products. Fixed physical sensors can therefore reverse mesh motion through
reference-coordinate inversion without differentiating point-location
topology.
All four arbitrary-order triangle vector bases also expose reference-coordinate
products. Monomial directional derivatives are evaluated analytically and
then transformed by the stored moment inverse. Because the active reference
coordinate has dimension two, each VJP uses exactly two analytical directional
evaluations; this avoids a tape and remains independent of polynomial-space
dimension in its number of sweeps.
RT, BDM, and first- and second-kind Nedelec physical-point observation wrappers
compose these basis products with affine inversion, Piola geometry, and local
coefficients. Their JVPs can move the sensor, element, and state simultaneously.
Their VJPs return sensor-position, element-vertex, and state-coefficient
cotangents, with the direct Piola shape term and inverse-coordinate shape term
accumulated exactly once.

Complex tetrahedral Nedelec PML elements expose the same geometry products
for all vertex coordinates, the complex Cartesian stretch, and wave number.
They compose the FortSym-generated curl--curl PML coefficient JVP/VJP with
the analytical Piola and quadrature products. Reverse mode returns real
geometry and wave-number gradients together with complex stretch cotangents
under the library's real-complex inner-product convention.

Global tetrahedral PML assembly preserves the primal merged CSC pattern in
its JVP. Its VJP accepts value cotangents aligned with that pattern, reverses
the Nedelec orientation transforms, and accumulates element geometry
cotangents at shared mesh vertices. Elementwise complex stretch gradients
remain separate, while the global wave-number gradient is summed across all
cells. This is the mesh-level interface consumed by sparse state adjoints.

The real tetrahedral curl--curl-plus-mass CSC assembler provides the
corresponding mesh, curl-coefficient, and mass-coefficient products. Its
reverse pass handles both first-order edge orientation and arbitrary-order
edge/face/cell basis transforms, then accumulates cotangents from adjacent
elements at shared vertices. Thus FEM--DtN and PML volume operators expose
the same geometry-facing contract.

Arbitrary-order tetrahedral H1 stiffness-plus-mass elements expose analytical
products for all vertex coordinates and both scalar coefficients. They
compose FortNum's FortSym-generated determinant and inverse products with
physical-gradient quadrature. One reverse sweep returns the complete element
shape gradient without forming the coordinate-to-matrix Jacobian.
The corresponding global CSC products preserve the merged arbitrary-order
DOF pattern. Reverse assembly gathers each sparse value cotangent directly
from CSC storage and accumulates element shape gradients at shared mesh
vertices, avoiding a dense global cotangent matrix.
The constrained tetrahedral H1 state API composes these products with fixed
Dirichlet elimination and an implicit sparse solve. Its JVP and VJP cover all
mesh coordinates, stiffness and mass coefficients, assembled volume loads,
and prescribed boundary values. The boundary mask remains an explicit
inactive topology input, while one reduced adjoint solve returns every active
cotangent needed by Poisson, diffusion-reaction, and scalar Helmholtz
optimization loops.

Scalar tetrahedral loads also have an explicit sampled-source contract.
Source values live at physical quadrature points and supplied spatial
gradients account for Eulerian point motion. A JVP combines mesh motion with
an independent source-parameter tangent; the VJP returns both shared vertex
cotangents and source-sample cotangents. This makes source parameters
composable without attempting to differentiate an opaque Fortran procedure
callback.
The sampled constrained-state wrapper composes this load path with global H1
assembly and fixed-mask elimination. One call differentiates an Eulerian
source, material coefficients, prescribed boundary values, and every mesh
coordinate through the complete Poisson or diffusion-reaction solve.
Fixed-reference tetrahedral H1 observations now expose the same product
contract. They differentiate the scalar value and physical gradient with
respect to every global solution coefficient and all vertices of the active
tetrahedron. The geometry path composes FortNum's FortSym-generated
three-by-three inverse products with the reference gradient, so a Poisson
state adjoint can continue through pointwise objectives without differentiating
an inverse or factorization numerically.
