---
title: Sheet-current representation parity
---

# Sheet-current representation parity

`evaluate_sheet_current_parity` is a small manufactured slab contract for
MHD-09.  It compares two representations of the same tangential surface
current without implementing an equilibrium or Maxwell solver:

* the explicit fitted/cut/DG/IGA surface ledger `A K`; and
* the resolved layer `int K delta_epsilon(d) dd dA` evaluated with the existing
  normalized Gaussian layer operator.

The contract requires one constant tangential vector `K` on the supplied normal
quadrature and reports both ledgers plus a scale-safe relative difference.  The
test uses an independently evaluated Gaussian oracle, so it is not a test that
merely reproduces the implementation.  Surface geometry, Ampere constitutive
laws, and PDE assembly remain caller-owned.

The slab fixture is the physical-first gallery seed for the later cylinder,
sphere, and torus parity cases.  It is intentionally one-dimensional in the
normal coordinate: geometry-specific surface parametrizations belong in
application repositories or higher-level examples.
