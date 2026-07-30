# Tetrahedral symbolic mixed Poisson

This example extends the symbolic Darcy form from triangles to a unit cube
split into six tetrahedra:

\[
 (q,\tau) - (p,\nabla\!\cdot\tau)=0,\qquad
 (\nabla\!\cdot q,v)=(1,v).
\]

It runs every implemented tetrahedral pair from RT0/DG0 through RT5/DG5.
The RT mass and rectangular divergence expressions compile directly into
FortSparse CSC blocks. For each tetrahedron, the constant DG test moment of
the numerical divergence must equal the independently computed geometric
volume. A second plot records the coupled system size as order increases.

CI generates `tetra_mixed_conservation_1d.png` and
`tetra_mixed_dofs_1d.png`; generated media are not committed.
