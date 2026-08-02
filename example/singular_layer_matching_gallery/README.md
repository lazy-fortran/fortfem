# Neutral singular-layer trace matching

This bounded MHD-08/MHD-09 foundation gallery samples analytical complex inner
and outer traces on a one-dimensional layer coordinate.  The first figure
shows their real profiles and the matched jump prescribed to
`assemble_singular_layer_matching`; complex components are retained in the
CSV source.  A second figure deliberately perturbs that jump and reports the
mismatch residual together with the exact JVP's centered-difference error.

The trace rows, states, jump, and positive weights are caller-owned data.  No
plasma closure, singular asymptotic, equilibrium reader, or production tearing
physics is selected.  Outputs are written to
`output/example/singular_layer_matching_gallery/`:

- `singular_layer_solution_1d.png`: physical outer/inner traces and matched jump;
- `singular_layer_diagnostics_1d.png`: mismatch and derivative diagnostics;
- `singular_layer.csv`: reproducible complex profile and residual samples;
- `provenance.json`: contract, analytical parameters, and diagnostic norms.
