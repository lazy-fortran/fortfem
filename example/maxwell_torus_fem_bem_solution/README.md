# Toroidal Maxwell FEM--BEM solution

This fixture solves the neutral curl--curl FEM--BEM system on a small exact
solid-torus mesh. A constant analytical (H(\mathrm{curl})) field supplies the
manufactured right-hand side, so the coupled solve has a known volume field
and zero surface-current solution. The first plot is the computed physical
field magnitude on a toroidal cross-section; a vector-arrow view, curved
geometry, error, and timing diagnostics follow.

The example exercises the same volume Nédélec, curved RWG trace, and dense
FEM--BEM system used by external Maxwell clients. It is a foundation fixture,
not a plasma or scattering model.
