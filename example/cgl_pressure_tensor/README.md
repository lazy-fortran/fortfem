# CGL pressure tensor example

This manufactured profile evaluates the gyrotropic CGL pressure tensor

\[
\mathbf P=p_\perp I+(p_\parallel-p_\perp)\mathbf b\mathbf b^T
\]

for a magnetic direction that rotates across a one-dimensional coordinate.
It also evaluates the generated product-rule force \(\nabla\cdot\mathbf P\),
checks the tensor trace independently, and writes reproducible CSV data.

The generated gallery plots show the diagonal/off-diagonal pressure
components and the two nonzero force-divergence components. This is a
constitutive manufactured-solution example, not a complete equilibrium or
transport solve.

Outputs:

- `cgl_pressure_components_1d.png`
- `cgl_force_divergence_1d.png`
- `cgl_profile.csv`
