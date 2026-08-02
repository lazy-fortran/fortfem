---
title: Toroidal spectral metadata contract
---

# Toroidal spectral metadata contract

`analyze_toroidal_spectral_modes` is the fixed-topology policy boundary for a
NESTOR-like modal trace.  Given a supplied list of `(degree, order)` labels and
a caller-selected rectangular window, it reports retained and omitted counts,
the explicit `(0,0)` zero-mode count, the coefficient energy of the supplied
list, and the energy of the modes outside that window.

The routine does not claim to estimate an unprovided spectral tail.  It also
does not silently remove a mean: `allow_zero_mode=.false.` returns a rejection
status when the supplied list contains `(0,0)`, leaving the zero-mode policy
visible to the caller.  This is the correct separation for free-boundary
solvers, where the mean constraint depends on the Green-kernel normalization
and on the chosen Neumann compatibility condition.

The energy diagnostics have analytical JVP/VJP actions for modal coefficients
with the real-part complex pairing.  The mode mask and counts are discrete
fixed-topology metadata and are therefore not differentiated.  The independent
test covers count/energy values, zero-mode rejection, central reassembly, and
the complex adjoint identity.
