# Generated mathematical kernels

The ordinary FortFEM library consumes checked-in Fortran and needs no CAS.
Mathematical definitions and symbolic derivatives live in `tools/codegen/app`;
`tools/codegen/check_generated.sh` regenerates and compares runtime outputs
against a clean FortSym revision declared in `tools/codegen/fortsym.lock`.

Generation prevents separately copied derivatives from drifting. It does not
prove the mathematical definition, element convention, shape validation,
assembly, quadrature, or floating-point behavior correct. Behavioral tests
therefore use independent interpolation, polynomial, conservation, and adjoint
oracles.

| Definition family | State | Independent evidence |
| --- | --- | --- |
| Triangle P1/P2 scalar reference values, gradients, Hessians | Generated | Exact polynomial reproduction, every P2 Hessian |
| Quadrilateral Q1 and interval P1 reference values/derivatives | Generated | Partition, polynomial reproduction |
| Triangle Whitney/RT0 values, curl/divergence | Generated | Oriented line/flux moments |
| Tetrahedron Whitney first-order values and curls | Generated | Oriented edge moments, rigid-rotation reproduction, all face Stokes fluxes |
| Arbitrary-degree triangle/tetrahedron scalar Lagrange jets | Generated | Degree 0–7 polynomial reproduction and signed AD products |
| Arbitrary-degree triangle first/second-kind Nedelec candidate jets | Generated | Exact oriented moments, degree 0–12 analytic monomial jets, zero-axis/vertex tangent and adjoint tests; RT/BDM inherit by rotation |
| Arbitrary-degree tetrahedron Nedelec/RT modal jets | Mixed generated/manual | Existing generated modal primitives; remaining derivative composition requires audit |
| Triangle Piola values, tangent and reverse products | Generated | Oriented line/normal moments, differentiated conservation, FD and adjoint products |
| Affine triangle/tetrahedron geometry and inverse-map products | Generated | Prescribed barycentric coordinates, joint-motion and translation invariance, large exact translations, FD and adjoint identities |
| Tetrahedron Piola values, tangent and reverse products | Generated | Oriented line and face flux moments, differentiated conservation, FD and adjoint products, invalid-geometry status regression |
| B-spline recurrence, polar/multipatch products | Mixed generated/manual | Geometry products generated; recurrence/operator audit pending |
| Nested mapped geometry, cut-cell moments and differential jets | Mixed generated/manual | Source-level inventory and focused oracle review pending |
| Weak operators and element/assembly products | Mixed generated/manual | Individual generated products exist; whole-owner audit pending |

The owner inventory at the pinned `91a7eb7` baseline had 475 runtime Fortran
files, including 60 generated files, and 52 generator programs. These counts are
an inventory, not a claim that every nongenerated file contains symbolic math.
Mesh topology, variable-degree iteration, shape/status checks, sparse assembly,
and linear algebra delegated to FortNum remain ordinary Fortran orchestration.

CMake currently exposes triangle compatibility and a C API subset. The scalar
arbitrary-order and tetrahedron Whitney interfaces are FPM-exposed; their modules
are not dependencies of those CMake targets. Add generated modules to any
explicit CMake source list when an exposed consumer begins using them.
