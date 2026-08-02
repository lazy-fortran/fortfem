# Field-aligned anisotropic 3-D diffusion

This neutral manufactured-solution gallery exercises the public three-dimensional
tensor-diffusion contraction with a strongly anisotropic field-aligned tensor.
The unit direction is a fixed oblique vector, so the large coefficient acts
along a direction that is not aligned with any coordinate axis:

\[
K = k_\perp I + (k_\parallel-k_\perp)\,b b^T,
\qquad b=(1,2,2)/3,
\qquad -\nabla\cdot(K\nabla u)=f.
\]

The exact solution is `u = sin(pi*x) sin(pi*y) sin(pi*z)` on the unit cube
with homogeneous Dirichlet data.  The tetrahedral P1 assembly evaluates `K`
through `evaluate_field_aligned_constitutive_tensor` and calls
`assemble_tensor_diffusion_matrix_3d`; no plasma model, closure, geometry
reader, or external data is involved.

The first figure is the computed 3-D solution with flux-direction arrows and
tetrahedral edges.  A subsequent mid-plane slice shows the scalar solution and
the in-plane anisotropic flux.  The executable writes an independent analytic
source/error/energy oracle to `diagnostics.csv`, plus a CSV nodal solution for
downstream plotting.

The constitutive construction follows the standard field-aligned projector
used in neutral anisotropic transport formulations.  The gallery is an
implementation fixture for future field-aligned FEM/IGA and not a plasma
equilibrium calculation.
