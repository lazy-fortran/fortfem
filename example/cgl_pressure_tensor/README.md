# CGL pressure tensor example

This manufactured profile evaluates the gyrotropic CGL pressure tensor

\[
\mathbf P=p_\perp I+(p_\parallel-p_\perp)\mathbf b\mathbf b^T
\]

for a magnetic direction that rotates across a two-dimensional sampled field.
It also evaluates the generated product-rule force \(\nabla\cdot\mathbf P\),
checks the tensor against an independent projector oracle, and writes
reproducible CSV data.

The physical-first gallery plot shows the tensor anisotropy and principal
direction field. The original one-dimensional diagonal/off-diagonal pressure
components and the two nonzero force-divergence components remain as
diagnostics. This is a constitutive manufactured-solution example, not a
complete equilibrium or transport solve.

Outputs:

- `cgl_tensor_principal_2d.png`
- `cgl_pressure_components_1d.png`
- `cgl_force_divergence_1d.png`
- `cgl_profile.csv`
- `cgl_tensor_field_2d.csv`
