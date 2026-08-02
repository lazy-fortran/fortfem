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

`test_force_balance_residual` compares the value against an independent
three-term contraction, checks a central-difference JVP and the adjoint
identity, and rejects a nonpositive volume measure.
