# Vector Maxwell FEM--BEM/DtN/PML parity

This gallery is a small, license-safe manufactured curl--curl fixture.  It
uses the same nonzero vector target to exercise three open-boundary paths:

* a curved toroidal Nédélec FEM--BEM coupling;
* the curved RWG/RBC Maxwell DtN map;
* Nédélec PML solves on a base box and a larger box.

The first images are physical fields and vector arrows.  Diagnostics (the
weak DtN response and the larger-domain target comparison) are written only
after the physical previews.  The benchmark is intentionally a foundation
fixture, not a plasma model and not a reproduction of external paper data.

The manufactured target is nonzero in all three components and includes a
constant rotational part.  The curved FEM--BEM solve reconstructs its edge
field on a toroidal cross-section; the DtN path applies the same trace-to-flux
map to a nonzero complex trace; and the PML path compares a base box with a
larger box.  `torus_vector_field.csv` and `pml_vector_field.csv` retain the
sampled physical vector values used by the plots.  `provenance.json` records
that the fixture is generated data rather than an external benchmark array.

For the independent fast gallery gate, set
`MAXWELL_VECTOR_OPEN_BOUNDARY_FAST=1`; this keeps the same operators and
manufactured target while using the smallest torus mesh.

Run it with:

```text
fo run --example maxwell_vector_open_boundary_parity_gallery
```

Outputs are written below
`output/example/maxwell_vector_open_boundary_parity_gallery/`.
