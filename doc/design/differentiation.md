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

Arbitrary-order tetrahedral Nedelec curl--curl-plus-mass element matrices
compose those Piola products with quadrature. Their JVPs cover all twelve
vertex coordinates and both material coefficients. A single analytical
reverse quadrature sweep returns the corresponding vertex and coefficient
cotangents, rather than requiring one forward sweep per design coordinate.
This is the preferred element path for large moving-mesh design spaces.

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
