# FCI interpolation-map construction

`build_fci_linear_interpolation_map_1d` is the first geometry-side map builder
for the PARALLAX-aligned FCI path.  Given strictly increasing source
coordinates and in-range target coordinates, it returns one row per target
with the two piecewise-linear weights that bracket it.  Endpoint rows are
represented exactly.

For source values (u_i), the map (M) satisfies the two basic independent
oracles

\[
  M\,1=1, \qquad M\,(a x+b)=a\,x_{\rm target}+b.
\]

The routine rejects non-finite or repeated source coordinates, targets outside
the closed source interval, and incompatible output shapes.  It is purposely
one-dimensional: a field-line service can use it for a coordinate slice, while
a future 2D/3D stencil builder can be verified against the same partition and
affine-reproduction contracts.  It does not copy PARALLAX geometry or lookup
code.

The focused test checks the exact affine oracle, partition of unity, endpoint
handling, and invalid-input paths.  Higher-dimensional interpolation,
geometry/map JVPs, and support-volume construction remain separate roadmap
items.

## Provenance

The separation between field-line tracing, interpolation maps, and support
operators follows the FCI architecture described by
[Stegmeir et al.](https://doi.org/10.1016/j.cpc.2016.12.014) and the local
PARALLAX documentation.  This implementation is an independent Fortran
contract; no PARALLAX source or benchmark data is included.
