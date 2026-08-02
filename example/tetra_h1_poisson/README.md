# Arbitrary-order tetrahedral H1 Poisson solve

This example solves

\[
  -\Delta u=f,\qquad
  u=0\ \text{on }\partial\Omega
\]

on a tetrahedron split into four conforming sub-tetrahedra. FortSym generates
the exact quartic bubble

\[
  u(x,y,z)=xyz(1-x-y-z)
\]

together with its gradient and source in one shared kernel. The public
FortFEM solver assembles degrees one through four directly into FortSparse
CSC storage, eliminates every boundary entity degree of freedom, and
reconstructs physical values and gradients.

The example writes uncommitted `convergence.csv` and
`tetra_h1_poisson_convergence.png` artifacts. Both \(L^2\) and gradient errors
must decrease under p-refinement, and the degree-four solution must reproduce
the quartic analytical field to roundoff. The physical view uses dense
tetrahedral samples with the conforming sub-tetrahedron edges overlaid, so the
solution is readable as a volume/slice view rather than an unlabelled cloud.
A small `gallery_sequence.txt`
record is also emitted as an execution oracle: it records the physical
solution plot before the convergence diagnostics are written.
