# Symbolic mixed Poisson

This example is the first vector finite-element problem after the scalar
Poisson introduction. It writes the mixed Darcy form directly:

\[
 (q,\tau) - (p,\nabla\!\cdot\tau)=0,\qquad
 (\nabla\!\cdot q,v)=(f,v),
\]

with \(q=-\nabla p\), RT1 fluxes, and discontinuous linear pressure. The
source and analytical solution are

\[
 p=\sin(\pi x)\sin(\pi y),\qquad
 f=2\pi^2\sin(\pi x)\sin(\pi y).
\]

The rectangular divergence blocks and their transpose are compiled from the
symbolic expressions into FortSparse CSC matrices. Three refinements verify
the expected second-order convergence of both flux and pressure. The first
plot is sampled from the solved DG1 pressure coefficients on the finest mesh;
the companion CSV records numerical and analytical values independently.

CI generates `mixed_poisson_solution_2d.png` first, followed by
`mixed_poisson_convergence_1d.png`; generated media are not committed.
