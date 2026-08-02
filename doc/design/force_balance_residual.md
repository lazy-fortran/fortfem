# Closure-neutral force-balance residual

`assemble_force_balance_residual` is the weak-composition boundary for
equilibrium and force-balance clients. For every test function and component
it adds three caller-owned contributions:

\[
R_{ic} = \sum_q w^V_q V_{qic}F^V_{qc}
       + \sum_b w^B_b B_{bic}F^B_{bc}
       + \sum_s w^S_s S_{sic}F^S_{sc}.
\]

The volume force can be supplied by magnetic, tensor-pressure, inertial, or
body-force blocks. Boundary traction and sheet-current terms are separate
inputs, so a distributional current is not silently regularized into a volume
load. Positive quadrature measures, finite data, and compatible shapes are
validated; zero-row boundary or sheet terms are allowed.

The JVP differentiates all test, force, and measure factors. The VJP is the
real dot-product adjoint of that product rule. Constitutive laws, coordinate
maps, units, pressure closures, Maxwell stress, Ampere jump laws, and physical
equilibrium variables remain caller-owned. This lets VMEC/GVEC/DESC-like,
SPEC-like, free-boundary, elasticity, and reduced-MHD clients share the same
weak residual without making FortFEM a plasma-physics code.

## Generated local contraction

The scalar product

\[
    c = w\,v\,f
\]

is emitted by `tools/codegen/app/gen_force_balance_products.f90` with FortSym
revision `fortsym@a0e3857436f15347a0ab69740e309d6598c64960`.  The committed
`src/generated/fortfem_force_balance_products.f90` contains value, JVP, and
VJP kernels; `force_balance_residual` calls the public wrapper
`evaluate_force_balance_product*` for volume, boundary, and sheet terms.  The
generator is included in `tools/codegen/generate.sh` and byte-checked by
`tools/codegen/check_generated.sh`, so a consumer build never needs a CAS and
the three contraction paths cannot drift independently.

`test_force_balance_residual` compares the value against an independent
three-term contraction, checks a central-difference JVP and the adjoint
identity, and rejects a nonpositive volume measure.

`test_force_balance_generated_product` checks the generated local product
against an independently written scalar oracle, a central difference, and the
real dot-product identity.  It intentionally does not parse generated source
or compare repository state.

## Shape-validation boundary

The volume block fixes the output shape: `volume_test` has shape
`(n_volume, n_test, n_component)`, `volume_force` has shape
`(n_volume, n_component)`, and `volume_weights` has length `n_volume`.
Boundary and sheet blocks use the same trailing dimensions and may have zero
rows when a caller has no such contribution.  Every measure must be finite and
strictly positive; an empty optional block has no measures and is accepted.
JVP increments must have exactly the corresponding primal shapes, while VJP
cotangents must have the corresponding output shapes.  A shape error returns
`FORTSPARSE_INVALID_MATRIX` before any contraction is attempted.

`test_force_balance_validation` exercises this boundary independently: it
checks the zero-row case, force and measure length mismatches, malformed JVP
increments, malformed VJP residual cotangents, and a compatible VJP output
allocation.  This keeps the neutral composition reusable for vector-valued
anisotropic, Maxwell, elastic, and sheet-current clients without relying on
Fortran's implicit shape assumptions.
