# Magnetic-axis isogeometric FEEC

This example constructs the Type-1 polar spline de Rham sequence

\[
H^1 \xrightarrow{\nabla} H(\mathrm{curl})
\xrightarrow{\mathrm{curl}} L^2
\]

on a periodic radial patch containing the magnetic axis. It assembles the
physical Piola-mapped Hodge operators and checks, for deterministic random
fields, that the discrete gradient and curl energies commute with their
corresponding Hodge maps.

The generated gallery artifacts lead with the computed scalar FEEC solution,
its physical gradient vectors, and the physical spline element lines. The
solution samples are dense in the periodic direction and include the magnetic
axis and outer boundary, so the first image shows the complete mapped patch
instead of a sparse parameter-space point cloud. A separate plot shows the same
curvilinear physical mesh without the field overlay. The energy identity errors
and measured assembly time are also reported. No rendered image is committed.

## Provenance

- A. Toshniwal and T. J. R. Hughes, *Isogeometric discrete differential forms:
  Non-uniform degrees, Bézier extraction, polar splines and flows on surfaces*,
  CMAME 376 (2021), 113576,
  <https://doi.org/10.1016/j.cma.2020.113576>.
- F. Holderied, S. Possanner, and X. Wang, *MHD-kinetic hybrid code based on
  structure-preserving finite elements with particles-in-cell*, JCP 433
  (2021), 110143, <https://doi.org/10.1016/j.jcp.2021.110143>.
