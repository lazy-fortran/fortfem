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
