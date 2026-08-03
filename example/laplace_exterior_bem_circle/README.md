# Exterior Laplace BEM on a circle

This example solves the unbounded exterior Dirichlet problem on the unit
circle with the manufactured decaying solution

\[
u(r,\theta)=\frac{\cos\theta}{r}.
\]

The periodic piecewise-constant single-layer density is assembled and solved
with FortFEM's dense solver, then evaluated directly at exterior points. The
gallery leads with the computed two-dimensional exterior solution, followed by
an analytical radial trace and the recovered boundary density. The interior of
the circle is shown as a neutral zero-valued region in the solution plot; no
artificial far boundary is introduced.

The run also writes `solution.csv`, `provenance.json`, and
`gallery_sequence.txt`.  The solution CSV contains the sampled computed field
and the exact value at each exterior point; the physical-first gate recomputes
the (x/r^2) oracle independently and checks the plot/data ordering.  Images
and CSV output are generated artifacts and are not committed.

## Provenance

- S. Sauter and C. Schwab, *Boundary Element Methods*, Springer, 2011,
  <https://doi.org/10.1007/978-3-540-68093-2>.
- FortFEM's logarithmic single-layer kernel uses the standard two-dimensional
  Laplace fundamental solution \(G(x,y)=-(2\pi)^{-1}\log|x-y|\).
