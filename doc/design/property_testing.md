---
title: Deterministic property testing
---

# Deterministic property testing

The `check` test helper now provides a small deterministic property path that
works with the ordinary `fo test` discovery flow. `property_rng_t` is a
per-case xorshift generator initialized from an explicit 32-bit seed; it does
not mutate Fortran's process-global `random_seed` state. Uniform and bounded
integer samples are reproducible across repeated runs with the same seed.

`check_property` generates an independent case seed for each requested case,
calls a callback that receives that seed, and reports the original seed when a
case fails. It repeatedly applies `property_shrink_integer` toward zero to
provide a compact deterministic witness. The callback can initialize its own
`property_rng_t` from the case seed, so rerunning one failure does not depend
on the order or number of random calls in previous cases.

The helper is intentionally small and test-only: it does not impose a
property framework or global random state on library clients. The independent
test checks reproducibility, bounds, shrinking, and the seeded runner through
the normal `fo` test path.
