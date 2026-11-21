---
title: Topology-adapted multipatch spline spaces
---

# Topology-adapted multipatch spline spaces

## Question

Which spline discretisation represents an equilibrium map and its fields on
poloidal planes that contain several magnetic axes, island chains, X-points
and a wall? It must have no polar singularity, a regularity chosen per knot,
block edge and vertex, local anisotropic refinement, and a fast matrix-free
energy Hessian. MILCA's frozen-map and non-nested lanes consume it in place
of C0 Lagrange P3 triangles.

## Model

A plane carries a layout of quadrilateral blocks. Blocks surround every
O-point (an O-grid without a pole). X-points are block vertices, and
separatrices, walls and the plasma boundary are block edges. Block k has a
ruled reference map

$$
G_k:\hat Q=[0,1]^2\to\hat\Omega,\qquad
G_k(\xi)=(1-\xi_2)\,C_k(\xi_1)+\xi_2\,L_k(\xi_1),
\tag{1}
$$

into the reference plane $\hat\Omega\subset\mathbb R^2$ with coordinates
$s=(s_1,s_2)$. $L_k$ is the straight segment between two corner vertices.
$C_k$ is either straight, which makes $G_k$ bilinear, or a boundary curve in
the block spline space. Interior interfaces are therefore always straight.
$G_k$ has control points $C_{k,i}+\gamma_j(L_{k,i}-C_{k,i})$ with the Greville
abscissae $\gamma_j$, which is exact by linear precision. On each block, $N^k_{ij}(\xi)=N_i(\xi_1;\Xi^1_k)
N_j(\xi_2;\Xi^2_k)$ are tensor B-splines of degree p on open knot vectors
$\Xi^d_k$. A breakpoint of multiplicity m gives $C^{p-m}$ there. Every global
basis function is a sparse combination of patch functions,

$$
\varphi_I\circ G_k=\sum_{ij}E_{I,(k,ij)}N^k_{ij}.
\tag{2}
$$

The space is fixed by the extraction matrix E:

| coupling | E |
|---|---|
| C0 across block edges | signed injection from `build_multipatch_signed_dof_map` |
| C1 across edges and at vertices | analysis-suitable G1 edge and vertex functions |
| THB refinement | truncated hierarchical functions in terms of fine-level patch functions |

Regularity is defined with respect to s. Along an interface the gluing data of
(1) depend only on the straight edge and on the transversal derivative, which
is linear along it, as for bilinear blocks. So $\alpha^{(1)}$, $\alpha^{(2)}$
and $\beta$ in (12) are linear, and the analysis-suitable G1 split of $\beta$
is checked per interface (see the C1 construction below). Curved
boundary edges are required: a C1 map from a polygon onto a smooth wall is
singular at the polygon corners, since its Jacobian cannot open a 90 degree
corner into a straight boundary. The remaining shape of walls, separatrices
and flux surfaces belongs to the unknown map $X(s)$, which lies in the same
space.

A field u with coefficients $c_{fI}$, $f=1..n_f$, in a toroidal basis
$T_m(\phi)$ (real Fourier modes or periodic B-splines along the curve chart) is

$$
u_f(s,\phi)=\sum_{m}\sum_{I}c_{f m I}\,T_m(\phi)\,\varphi_I(s).
\tag{3}
$$

The energy is a sum of pointwise kernels over cells e, quadrature points q and
toroidal planes j,

$$
\mathcal E(c)=\sum_{e,q,j}w_{qj}\,W_{eqj}\bigl(J_{eqj}(c)\bigr),
\qquad J_{eqj}(c)=\Phi_{eqj}c,
\tag{4}
$$

where the jet $J$ holds the values and derivatives of all fields up to
order 2 at the point. Component $f+n_f d$ is field f with d = 0 for the value
and d = 1, 2, 3 for $\partial_{s_1},\partial_{s_2},\partial_\phi$, followed by
the second derivatives in the order (11, 21, 31, 22, 32, 33). The weight $w_{qj}$
includes $\det DG_k$ and the plane spacing.

## Derivation

Cell evaluation. On element $e=(e_1,e_2)$ of block k with breakpoint widths
$h_1,h_2$ and local variable $t\in[0,1]^2$, the 1D Bezier extraction (Borden,
Scott, Evans, Hughes 2011, Algorithm 1) gives

$$
N_{f_{e_d}+a-1}(\xi_d)=\sum_{b=0}^{p}C^d_{e_d}(a,b)\,B^p_b(t_d),
\tag{5}
$$

$$
\varphi_I\big|_e=\sum_{a,b}\bigl(E\,(C^1_{e_1}\otimes C^2_{e_2})\bigr)_{I,ab}
\,B^p_a(t_1)B^p_b(t_2),
\tag{6}
$$

where $B^p_b(t)=\binom pb t^b(1-t)^{p-b}$ and its first two derivatives are
FortSym kernels (`gen_interval_bernstein_jets`). One evaluation path serves
every regularity, because only E changes.

Derivatives with respect to s. With $\xi=\xi_0+h\odot t$ and $s=G_k(\xi)$,

$$
\nabla_s u=DG_k^{-T}\,\nabla_\xi u,
\tag{7}
$$

$$
\nabla_s^2u=DG_k^{-T}\Bigl(\nabla_\xi^2u-\sum_{r}(\nabla_s u)_r\,
\nabla_\xi^2G_{k,r}\Bigr)DG_k^{-1}.
\tag{8}
$$

$DG_k$ and $\nabla_\xi^2G_k$ come from the same cell evaluation as u, applied
to the control points of (1).

Gradient and Hessian. J is linear in c, so

$$
\partial_c\mathcal E=\sum w\,\Phi^T\partial_JW,
\tag{9}
$$

$$
\partial_c^2\mathcal E=\sum w\,\Phi^T\,\partial_J^2W\,\Phi.
\tag{10}
$$

The kernel returns $\partial_J^2W$ dense or packed (LAPACK 'U' column-major).
For Fourier modes the sum over planes in (10) couples modes m, m' through
$T_m(\phi_j)T_{m'}(\phi_j)$ and their derivatives. The library contracts it,
and the kernel only sees points (cell, plane).

Sum factorisation. On a cell the Bernstein coefficients are
$u^e=(E\,C_e)^Tc$. Values and derivatives at $(p+1)^2$ Gauss points follow from
two 1D contractions with the Bernstein tables, at a cost of $O(p^3)$ per cell
and field instead of $O(p^4)$. The Hessian action applies (10) as
gather, two contractions, the pointwise product, two transposed contractions
and a scatter, without forming the matrix. The assembled CSR Hessian is the
reference.

Knot grading. Breakpoints come from the GVEC radial grid types (`sgrid.F90`,
MIT): uniform, $\sqrt s$ (finer at s = 1), $s^2$, bump, and bump-edge.
Explicit breakpoints with multiplicities lower the regularity exactly at
physical kinks.

## Gates

Each gate runs in seconds. CPU-heavy ladders run on sCluster.

1. C1 space on the 5-block O-grid disc (centre square, four trapezoids, four
   valence-3 vertices at mid-radius). Poisson with an exact solution
   converges at order p + 1 in L2 for p = 3, 4, 5. Jumps of the first
   derivatives across block edges and at the vertices are at round-off.
2. THB refinement towards a boundary layer at the wall reaches a given error
   with fewer DOFs than uniform knots. Target: a cell width of 0.005 in s at
   the wall (DIII-D edge layer).
3. Two-O-point layout (doublet-like) and a layout with an X-point vertex:
   mesh and solve.
4. Curve-chart 3D extension on a figure-8 guide curve.
5. De Rham complex on the multipatch layouts: cohomology dimensions equal the
   Betti numbers.
6. Export of piecewise-linear A on tetrahedra (lowest-order Nedelec edge
   integrals, so div B = 0 per tetrahedron exactly) and of smooth tensor data.
   The projection error converges.

## Delivery

M1: graded knots, Bezier extraction, Bernstein jets, explicit layouts and the
5-block O-grid, C0 space, jets to order 2, energy hooks, O-grid Poisson.
M2: C1 (gate 1). M3: THB (gate 2). M4: layout generator (gate 3).
M5: toroidal direction (gate 4). M6: de Rham (gate 5). M7: exports (gate 6).

## C1 construction

The C1 space is the Argyris-type space of Kapl, Sangalli and Takacs
(arXiv:1711.05161, CAGD 2019, and the survey arXiv:1812.09088) on the
reference blocks (1). It is a direct sum of patch-interior, edge and vertex
functions, and its dimension does not depend on the geometry. For degree p,
regularity r, n elements per direction and $N=(p-r)(n-1)+p+1$ functions per
direction,

$$
\dim A=|I_\Omega|\,(N-4)^2+|I_\Sigma|\,\bigl(2(p-r-1)(n-1)+2p-9\bigr)+6|I_\chi|,
\tag{11}
$$

where $I_\Omega$, $I_\Sigma$ and $I_\chi$ are the blocks, interior edges and
vertices. The edge term is the count over the kept edge indices. Both papers
print $p-9$, which contradicts that count and the Bezier case p = 5, n = 1
(one edge function). Each vertex of any valence carries six functions that
interpolate $(\varphi,\nabla\varphi,\nabla^2\varphi)$ there.

On an interface in standard form, $G^{(1)}(0,\xi)=G^{(2)}(\xi,0)$, the gluing
data on the edge are

$$
\alpha^{(1)}=\det[\partial_1G^{(1)},\partial_2G^{(1)}],\quad
\alpha^{(2)}=\det[\partial_1G^{(2)},\partial_2G^{(2)}],\quad
\beta=\det[\partial_2G^{(2)},\partial_1G^{(1)}],
\tag{12}
$$

with $\beta=\alpha^{(1)}\beta^{(2)}+\alpha^{(2)}\beta^{(1)}$ and linear
$\beta^{(i)}$. With $b^+_j\in S^{p,r+1}$, $b^-_j\in S^{p-1,r}$, $c_0=b_0+b_1$
and $c_1=(h/p)b_1$, the edge functions on block 1 are

$$
f_{j,0}=b^+_j(\xi_2)c_0(\xi_1)-\beta^{(1)}(\xi_2)\,{b^+_j}'(\xi_2)\,c_1(\xi_1),
\qquad
f_{j,1}=\alpha^{(1)}(\xi_2)\,b^-_j(\xi_2)\,b_1(\xi_1),
\tag{13}
$$

and on block 2 with $\xi_1\leftrightarrow\xi_2$ and the opposite sign of
$f_{j,1}$. Each function is a combination of patch B-splines and enters E in
(2) as one sparse row.

Constraints from the papers:

- $p\ge3$ and $1\le r\le p-2$. A C1-coupled layout therefore has interior knot
  multiplicity at least 2, while C0-coupled layouts keep $C^{p-1}$.
- AS-G1 is the only class with optimal h-approximation for $r\le p-2$
  (Collin, Sangalli, Takacs, CAGD 2016). Blocks with curved interior
  interfaces would lose order, so interior curvature is carried by the map.
  The linear $\alpha^{(i)}$, $\beta^{(i)}$ are fitted on each interface, and a
  nonzero fit residual rejects the layout.
- At r = p - 2 the basis is locally linearly dependent. THB refinement then
  relies on condition (P1) and the vertex-patch marking of Bracco, Giannelli,
  Kapl, Vazquez (arXiv:2204.10000).

THB-splines follow Buffa, Gantner, Giannelli, Praetorius, Vazquez
(arXiv:2107.02023): truncation, T-admissible marking of class mu and Dorfler
marking. On C0 layouts the estimator needs edge-jump terms. The C0 de Rham
complex uses broken FEEC with averaging conforming projections (Guclu,
Hadjout, Campos Pinto, arXiv:2208.05238; psydac
`feec/conforming_projectors.py`, MIT). Knot operations follow igakit
`igalib.f90` (BSD-2). G+Smo `gsUnstructuredSplines` and `gsHSplines` (MPL-2)
were read for their data structures only.

## Open

- Mixed regularity, $C^{p-1}$ inside blocks with $C^{p-2}$ only near C1
  interfaces, would keep the higher interior smoothness. The papers read
  assume uniform r.
- The vertex-function scaling of arXiv:1812.09088 (p. 13) is to be checked
  against the original before M2.
- Layout generation from O-points, X-points and curves.
