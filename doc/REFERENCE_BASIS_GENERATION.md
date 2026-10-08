# Reference basis kernels

The six legacy reference-element modules consume checked-in FortSym-generated
kernels. Normal library builds do not need FortSym or Python. The Fortran source
`tools/codegen/app/gen_reference_basis_products.f90` owns their mathematical
basis functions and derives gradients, Hessians, curl and divergence from them.
Geometry transforms, Piola maps and assembly remain separate operations.

| Family | Generated quantities |
| --- | --- |
| Triangle P1/P2 | Values, gradients and Hessians |
| Quadrilateral Q1 | Values, gradients and Hessians |
| Interval P1 | Values and derivatives |
| Triangle Whitney edge | Reference vectors, curl and divergence |
| Triangle RT0 | Reference vectors, divergence and curl |

The existing public APIs consume the quantities they expose; unused quantities
stay within the generated jet. The edge divergence routine describes the
reference vectors. A general physical covariant-Piola divergence needs the full
metric, and cannot be specified by triangle area alone.

The previous handwritten P2 Hessian omitted the `xi,xi` derivative of
`4*(1-xi-eta)*xi` and the `eta,eta` derivative of
`4*eta*(1-xi-eta)`. The former test checked only the first vertex basis function.
An independent polynomial-reproduction test fails 15 of 18 checks on the
parent: constant and affine Hessians are also corrupted by these omissions.
It now checks all constant, affine and quadratic polynomials at interior and
boundary points. Independent circulation/normal-flux and bilinear reproduction
checks cover the remaining migrated reference families.

The code-generation FortSym pin includes the owning SymEngine ABI repair.
The previous pin cannot compile against SymEngine 0.15's opaque C expressions.
This changes the offline generator only; it does not add a runtime dependency.
Generate with `tools/codegen/generate.sh`, and compare reproducible outputs with
`tools/codegen/check_generated.sh`. Generation is reproducibility evidence;
independent polynomial and integral checks establish behavior.

Arbitrary-order triangle/tetrahedron recurrences are separate from these six
legacy modules. Their derivative-product leaves need the generic generated-jet
migration; this document does not claim those paths have already been replaced.
