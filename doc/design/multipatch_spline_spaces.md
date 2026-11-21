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

Regularity is defined with respect to s. On bilinear blocks the gluing data
(12) are linear, so the blocks are analysis-suitable G1 and the C1 construction
below applies unchanged. Ruled blocks with an arc side are used for C0 spaces
only. For them the transversal derivative across the interface opposite the
arc depends on the arc, and at the inner corners the mixed derivative would
force the arc tangent onto its chord.

A C1 map cannot meet a smooth wall at a corner of a polygonal reference. If
every interface at a reference corner is C1, then DX is continuous there, and
X maps the corner's reference angle (below 180 degrees) onto a straight wall
only with det DX = 0. C1 layouts therefore use per-edge regularity.
Interfaces that end at a reference corner where a smooth piece meets the
polygon are C0, which lets DX jump at that corner. All other interfaces and
interior vertices are C1. The remaining shape of walls, separatrices and flux
surfaces belongs to the unknown map $X(s)$, which lies in the same space.

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

## Per-edge regularity

Interior sides flagged C0 are coupled in C0 inside a C1 space. On such a
side both blocks use the boundary gluing data ($\alpha=1$, $\beta=0$), the trace
functions $f_{j,0}$ of (13) are shared, and the transversal functions
$f_{j,1}$ belong to each side. At a vertex the corners joined by C1 sides form
sectors S. Each sector carries C2 data $\phi_S=(\varphi,\nabla\varphi,\nabla^2\varphi)$.
Across a C0 side with tangent t between sectors S and S' the traces agree
when

$$
\varphi_S=\varphi_{S'},\qquad \nabla\varphi_S\cdot t=\nabla\varphi_{S'}\cdot t,\qquad
t^T\nabla^2\varphi_S\,t=t^T\nabla^2\varphi_{S'}\,t,
\tag{18}
$$

since the sides of bilinear blocks are straight. The vertex functions are
the functions (13) built from an orthonormal basis of the null space of (18),
6 per vertex with one sector and $1+3c$ per vertex with c C0 sides
separating c sectors (c = 2 at a boundary vertex: 9).

A curved side forbids C1 on its two adjacent sides. A C1 space uses the map
$G=B+(1-\tau)^2(C-\text{chord})$ (B bilinear, $\tau$ the distance parameter from the
curved side C), whose transversal derivative on the opposite side equals
that of B. The interface opposite C therefore keeps the linear gluing data
of B and stays AS-G1. In the O-grid the ring-ring sides are C0 and the
centre-ring sides C1.

The C1 core meets the periodic annulus in C0. The core functions with a
trace on the interface (shared edge functions $f_{j,0}$ and vertex functions
at the interface vertices) are replaced by (i) an orthonormal basis of the
zero-trace combinations of the vertex functions at each interface vertex,
and (ii) for each annulus interface function, the minimal-norm local
combination of core functions whose trace equals its restriction (16). In a
C1 space the extraction depends on the vertices, so a geometry update
rebuilds it with the same dimension.

Tests (p = 3, u - lap u = f, zero normal flux on the unit disc):

| layout | n = 4 | n = 8 | n = 16 |
|---|---|---|---|
| curved O-grid, C1 with C0 ring-ring sides | 5.0e-4 | 3.4e-5 | 2.3e-6 |
| C1 core + periodic annulus | 1.2e-3 | 7.3e-5 | 4.5e-6 |

On the polygonal O-grid with C0 ring-ring sides the dimension is the
Kapl-Sangalli-Takacs count plus $n_1-4$ per C0 side and 3 per split
boundary vertex. E has full row rank, values are continuous, gradient
jumps across C1 sides are below 1e-11, and across C0 sides they are O(1).
After a new D-shaped interface the C1 core stays C1 (without the rebuild
the check fails).

## De Rham complex

On a block with B-splines $N_i$ of degree p and the derivative splines
$D_i=\frac{p}{t_{i+p+1}-t_{i+1}}N_{i+1,p-1}$ (Curry and Schoenberg), so that
$N_i'=D_{i-1}-D_i$, the section forms are

$$
N_iN_j,\qquad D_iN_j\,d\xi_1,\ N_iD_j\,d\xi_2,\qquad D_iD_j\,d\xi_1\wedge d\xi_2 .
\tag{22}
$$

Their degrees of freedom are the nodes, edges and faces of the control mesh.
Nodes are the C0 functions. Interface edges are shared, each oriented from the
lower to the higher global node with a sign per block, which makes the
tangential trace of 1-forms continuous. Faces are broken. grad and curl are
incidence matrices, independent of the geometry, and the edge graph is the one
used for tree-cotree gauges. Along t the nodal basis $L_k$ is the
trigonometric cardinal basis at $t_k=(k-1)/n$, $n=2M+1$, and the edge basis
$H_k$ the trigonometric histopolation basis on $[t_k,t_{k+1}]$, with
$L_k'=H_{k-1}-H_k$:

$$
H_k(t)=1+2\sum_{m=1}^{M}\frac{\cos 2\pi m(t-t_k-h/2)}{\operatorname{sinc}(m/n)},\qquad h=1/n .
\tag{23}
$$

The 3D spaces are $V_0=N\otimes T$, $V_1=[E\otimes T;N\otimes T']$,
$V_2=[F\otimes T;E\otimes T']$, $V_3=F\otimes T'$ with
grad $=[G\otimes I;I\otimes D_t]$, curl $=[[C\otimes I,0];[I\otimes D_t,-G\otimes I]]$,
div $=[I\otimes D_t,-C\otimes I]$. The degrees of freedom of the projections
are values at the Greville points, integrals along the Greville edges and
integrals over the Greville cells (times t-nodes or t-intervals). By Stokes
they satisfy $d\,\Pi_k=\Pi_{k+1}d$ up to the quadrature of these
integrals. Fields enter as parametric proxies, and the user applies the
covariant or Piola map with $J=\partial x/\partial(\xi,t)$.

Tests (gate 5):
- Betti numbers from ranks: disc (1, 0, 0), ring (1, 1, 0),
  disc $\times$ circle (1, 1, 0, 0), ring $\times$ circle (1, 2, 1, 0). d d = 0
  exactly.
- Commuting defects: 2D grad and curl below 1e-10 (p = 3, n = 4). 3D on a
  torus with the cylindrical chart (p = 2, M = 2): 1.4e-14 (grad), 4.3e-14
  (curl), 1.2e-14 (div).
- Projection errors (p = 3, n = 4 and 8): 0-forms 2.9e-5 and 3.6e-6,
  1-form proxies 3.0e-6 and 3.5e-7.

Gauge strategy (agreed with MILCA): ungauged Krylov solves with an
auxiliary-space (Hiptmair-Xu) preconditioner. The energy depends on A only
through curl A, so its gradient is orthogonal to ker(curl) and every
Newton system is consistent. Tree and cotree serve the cohomology, the
boundary data and gauged direct solves of small cases.

- Tree: a breadth-first spanning forest of the edge graph of grad.
- First cohomology: section classes come from a dual spanning tree of faces
  and cotree edges rooted at the outside. Every leftover cotree edge
  generates one closed form, completed leaf first along the dual tree. In
  3D the section classes are repeated on every t-node, and the toroidal
  class is one on all t-edges. All representatives are shifted by a
  gradient to vanish on the tree.
- Boundary surface (curve times circle): its 1- and 2-form dofs, the
  surface curl (restriction of curl, for gauge-invariant B·n coupling)
  and the toroidal and poloidal surface classes (the two flux numbers).
- Mass matrices of all forms over the parametric cells (reference
  measure), with a caller-supplied pullback weight for the
  parametric proxies.

Tests: tree size (vertices - 1), b1 representatives for disc, ring, disc
$\times$ circle and ring $\times$ circle, curl h = 0, h = 0 on the tree,
rank([G, h]) = rank(G) + b1, surface curl of gradient traces = 0, two
independent surface classes, symmetric mass matrices reproducing
independently integrated norms (2D and 3D).

The complex needs a C0 section whose extraction is a signed injection
(O-grid, X-point and curved layouts). Not yet covered: the annulus with
its interface constraint (the 3D lane uses the leg-free X-point layout,
a conforming 4-block annulus, instead) and THB sections.

Periodic splines along t: nodal basis $N_k$ of degree $q\ge2$ on $n$
uniform elements, edge basis $D_k=N_{k,q-1}/h$ with $N_k'=D_{k-1}-D_k$, so
$d_t$ is the same cyclic difference. The degrees of freedom are values at the
Greville points $t_k$ and integrals over $[t_k,t_k+h]$, and the projections
apply the inverse periodic collocation or histopolation matrix along t. Tests:
b1 representatives of disc and ring $\times$ circle, mass matrices,
$\Pi$ on constants, and 3D commuting defects of 2.5e-14 (grad),
6.4e-14 (curl), 1.9e-14 (div) with q = 3, n = 5.

The auxiliary-space interpolation $\Pi$ maps nodal vector splines
$\varphi_i e_d$ (Cartesian components) to the commuting projection of
their 1-forms, with proxies $J^Te_d\varphi_i$. Per block it uses the 1D
inverses of the collocation and histopolation matrices, truncated at a
relative tolerance (default 1e-6; the entries decay geometrically). Test:
$\Pi$ applied to a constant field $e_d$ equals $G\,\Pi_0x_d$ (2D and 3D).
Nonzeros per column on the torus (p = 2, M = 2, n = 2, 4, 8): 26, 56, 73 at
1e-6 and 26, 29, 26 at 1e-4. The growth at 1e-6 ends once a block is longer
than the truncation band.

## Toroidal direction

Fields are $u_f(s,t)=\sum_{m,I}c_{fmI}T_m(t)\varphi_I(s)$ (3) with t in [0, 1)
covering one field period. $T_m$ are the real Fourier modes $1,\cos 2\pi mt,
\sin 2\pi mt$ ($m\le M$) on $J\ge 2M+1$ equispaced planes, or periodic
splines of degree q with Gauss planes. A jet is formed on each plane j from the
section s-jets of every mode,

$$
\partial_s^\alpha\partial_t^k u_f(s_q,t_j)=\sum_m T_m^{(k)}(t_j)\,
\bigl(\partial_s^\alpha u_{fm}\bigr)(s_q),
\tag{21}
$$

and the kernel sees one (cell, plane) at a time with weight $w_qw_j$. The
gradient and the Hessian apply the transposed contraction. The assembled Hessian
is block sparse: section dof pairs that share a cell couple through dense
$(n_f n_m)^2$ blocks.

A chart maps (s1, s2, t) to x. Cylindrical:
$x=(s_1\cos\phi,s_1\sin\phi,s_2)$, $\phi=2\pi t/n_{fp}$. Guide curve:
$x=r(t)+s_1e_1(t)+s_2e_2(t)$, with r the trigonometric interpolant of the
samples and $(e_1,e_2)$ the rotation-minimising frame (double reflection on
4096 points), rotated back by its holonomy angle times t and kept to 64
harmonics. Kernels apply the metric $J=\partial x/\partial(s,t)$ themselves,
e.g. $W=(\tfrac12|J^{-T}\nabla u|^2+\tfrac12u^2-fu)|\det J|$.

Tests (p = 3, O-grid section, $u-\Delta u=f$, exact Neumann data):

| case | toroidal basis | L2 errors (n = 2, 4, 8) | rate |
|---|---|---|---|
| torus, $n_{fp}=1$, $u=x_1+x_2^2+x_1x_3$ | Fourier, M = 2 | 5.2e-4, 4.5e-5, 2.4e-6 | 4.2 |
| torus, $n_{fp}=2$, $u=x_1^2-x_2^2+x_3$ | Fourier, M = 2 | 2.0e-4, 2.7e-5, 1.4e-6 | 4.2 |

- Torus with periodic cubic splines in t on 8, 16 and 32 elements (n = 4):
  5.9e-2, 2.1e-3 and 1.2e-4 (rate 4.2).
- Figure-8 tube (gate 4): guide curve $(2\cos\theta,\sin2\theta,1.2\sin\theta)$,
  radius 0.3, $u=\sin(x_1/2)\cos(x_2/2)+0.3x_1x_3$, n = 2. With M = 5, 8 and
  11 the errors are 6.3e-3, 1.3e-3 and 4.8e-4. The frame of this curve has
  slowly decaying harmonics, so the Fourier error dominates at these M.
- The assembled Hessian equals the matrix-free product to 1e-15 relative,
  and the gradient matches central differences to 3e-9 relative.

## X-point layouts

The plasma boundary $\Gamma_p$ is given as four sampled segments between
the corners $K_1$ (bottom), $K_2$ (outboard), $K_3$ (top) and $K_4$ (inboard). A
corner with two sampled legs to the wall is an X-point. The layout consists of:

| region | blocks |
|---|---|
| core: O-grid, inner square $c+\beta(K_q-c)$ | 5 |
| scrape-off layer of each segment, inner and outer layer split at the mid points of its radial sides | 8 |
| private flux region of each X-point: three quads around $(m_a+m_b+m_w)/3$, with $m_a$, $m_b$ the leg mid points and $m_w$ the wall mid point | 3 per X-point |

The radial sides are the legs at X-points and the straight ray from c
through $K_q$ at smooth corners. Segments and legs are clamped splines
through their samples at chord length (knots by averaging, Piegl and Tiller,
eq. 9.8). The last leg sample is moved to the closest wall point. Wall
pieces are windows of the periodic wall spline. Every block map is the Coons
patch

$$
G(u,v)=(1-v)S_1(u)+vS_3(u)+(1-u)S_4(v)+uS_2(v)-B(u,v)
\tag{20}
$$

of its four sides $S_i$ (B bilinear in the corners), interpolated at the
Greville points. The boundary rows are the Greville interpolants of the
sides, so the maps of neighbouring blocks agree. All interfaces are C0. The core side of
$\Gamma_p$ carries LAYOUT_PLASMA_BOUNDARY, the scrape-off sides of the legs
LAYOUT_SEPARATRIX. The same corners with legs give the same topology, so a
new shape keeps the extraction. A geometry update to a different topology
returns status 1, and fields are then moved by projection.

Test (p = 3, double null, $\Gamma_p$ of quadratic Bezier arcs with 90 degree
corners at the X-points, elliptic wall), $u-\Delta u=f$ with exact Neumann
data:

| n | L2 error | area error | $\oint_{\Gamma_p}R\,dZ$ error |
|---|---|---|---|
| 2 | 3.1e-4 | 2.7e-4 | 3.6e-4 |
| 4 | 4.7e-5 | 7.6e-5 | 9.9e-5 |
| 8 | 2.5e-6 | 6.1e-6 | 8.6e-6 |
| 16 | 1.6e-7 | 3.9e-7 | 5.2e-7 |

A moved double null (other X-point positions, elongation and radial shift)
keeps the space and its accuracy. The single-null layout has 16 blocks, and
updating a double-null space to it returns status 1.

A secondary X-point inside the scrape-off layer (near double null) is given
at a smooth corner. It becomes the vertex where the two layers of that
corner split, so the C0 lines of the layout pass through it. Radial THB
levels refine the radial direction of the core ring and scrape-off blocks
and both directions of the private flux blocks, whose sides lie on the legs,
which keeps the levels conforming. On the double null with a bump of width
0.08 at the lower X-point (exact H1 error as indicator, two radial levels,
then isotropic ones), THB reaches an H1 error of 6.2e-2 with 2686 DOFs,
against 6.7e-2 with 6229 DOFs for uniform knots.

## THB refinement

Level spaces $V_0\subset V_1\subset\dots$ are the C0 multipatch spaces (with
the annulus constraint) on dyadically refined breakpoints: all blocks, or
only the radial knots of the annulus (anisotropic, chosen when a level is
created). Every level map is the exact knot insertion of the level-0 map.
The domains $\Omega_m$ are unions of level-m cells, with $\Omega_0$ the
whole plane section. With $\mathcal B_m$ the global basis of $V_m$, the
active functions and the truncation are

$$
\mathcal A_m=\{\varphi\in\mathcal B_m:\ \mathrm{supp}\,\varphi\subseteq\Omega_m,\
\mathrm{supp}\,\varphi\not\subseteq\Omega_{m+1}\},\qquad
\mathrm{trunc}^{m+1}\Bigl(\sum_{\psi\in\mathcal B_{m+1}}c_\psi\psi\Bigr)=
\sum_{\mathrm{supp}\,\psi\not\subseteq\Omega_{m+1}}c_\psi\psi ,
\tag{19}
$$

and the THB basis is $\{\mathrm{trunc}^{L}\cdots\mathrm{trunc}^{m+1}\varphi:\ \varphi\in\mathcal A_m\}$
(Giannelli, Jüttler, Speleers 2012). The coefficients of $\varphi\in\mathcal B_m$
in $\mathcal B_{m+1}$ come from tensor knot insertion of its broken functions, read
off at one broken function per fine function whose extraction column holds
that function alone. On the active level-n cells a truncated function equals its
level-n representation after the truncations up to n. These representations are
its extraction rows, so the hierarchical space is an ordinary spline space
whose patches are (block, level) pairs and whose cells are the active cells.
Marking is Dörfler marking of a caller-supplied indicator. Refinement adds
the class-2 admissible closure (Buffa, Gantner, Giannelli, Praetorius,
Vázquez, arXiv:2107.02023): refining a level-m cell Q also refines the active
level-(m-1) cells in the support extension of Q's parent. Every truncated
function then spans at most two levels of active cells. Fields move between
spaces by L2 projection, which is exact for nested spaces.

Test (p = 3, O-grid with periodic annulus, base n = 4):

- After three local refinements (the second radial), the functions with all
  coefficients 1 sum to 1 to 1e-12. A level-0 field and a coarser
  hierarchy field are reproduced by projection to 1e-10. No function
  spans more than two levels.
- Gate 2: wall layer $u=e^{-(1-r)/\epsilon}e^{\kappa(\cos(\phi-\phi_0)-1)}$,
  $\epsilon=0.04$, $\kappa=6$, Neumann data, exact cellwise H1 error as
  indicator, Dörfler $\theta=0.5$, two radial levels, then isotropic ones:

| space | DOFs | H1 error |
|---|---|---|
| uniform n = 8 | 833 | 2.5e-2 |
| uniform n = 16 | 2801 | 3.4e-3 |
| THB, step 6 | 602 | 1.3e-3 |

A refinement step whose marked region is smaller than one level-(m+1)
support adds cells but no functions (THB step 2: 281 DOFs unchanged).

## Periodic annulus

The single-axis layout is a periodic annulus between the wall and an
interface circle, around the five-block O-grid core. Along $\xi_1$ (the
angle, counter-clockwise from the core corner) the annulus carries the
$n+p$ unclamped B-splines $N_i$ on the periodically extended knots
$x_{i-n}-1,\;x_i,\;x_{i+n}+1$. The functions $i$ and $i+n$ coincide on
[0, 1], and E identifies them:

$$
\varphi_{(i,j)}\circ G_a=N_i(\xi_1)N_j(\xi_2)+[\,i\le p\,]\,N_{i+n}(\xi_1)N_j(\xi_2),
\qquad i=1..n.
\tag{15}
$$

The $\xi_1$ breakpoints are the angles of the breakpoints of the core sides
on the interface. The core side parameter t is affine in the angle, so the
restriction of $N_i$ to a core side lies in that side's spline space and

$$
N_i\bigl(\xi_1(t)\bigr)=\sum_m R_{mi}\,B_m(t)
\tag{16}
$$

is exact (Greville interpolation). The core functions $B_m$ on the
interface are replaced by $\sum_m R_{mi}B_m$ attached to the annulus trace
function i, so values are continuous (C0) and the core trace lies in the
annulus trace space. The annulus map is $c+\rho(\xi_2)\,r(\xi_1)$ with r
the periodic Greville interpolant of the unit circle and $\rho$ linear from
the wall ($\xi_2=0$) to the interface ($\xi_2=1$). The core sides on the
interface take $c+\rho_{\rm int}\sum_m R_{mi}r_i$ as their control points,
so both maps agree exactly. The radial knots of the annulus are free, e.g.
graded toward the wall.

Test (p = 3, wall radius 1, interface 0.6, core corners 0.3, n cells per
core block direction and annulus radius, 4n in angle), u - lap u = f with
zero normal derivative:

| n | L2 error | area error |
|---|---|---|
| 2 | 2.1e-2 | 3.8e-3 |
| 4 | 1.2e-3 | 2.2e-4 |
| 8 | 7.3e-5 | 1.3e-5 |

Seam jumps of value, gradient and Hessian are below 1e-9 relative. Interface
jumps of value and map are below 1e-12.

## Closed curves

Wall and interface may be closed curves sampled at points $S_q$, q = 1..N,
counter-clockwise from the outboard midplane. With the normalised chord
length $u_q=\sum_{l<q}|S_{l+1}-S_l|/L$ the curve is the periodic spline of odd
degree $p_c$ (default 3) with knots $u_q$ that interpolates the samples,

$$
C(u)=\sum_{i=1}^{N}P_iN_i(u),\qquad C(u_q)=S_q .
\tag{17}
$$

Every side carries a fixed parameter window: side point
$C(u_{\rm off}+u_0+\Delta u\,t)$, with $u_{\rm off}=5/8$, $u_0=0,\tfrac14,\tfrac12,\tfrac34$
and $\Delta u=\tfrac14$ for the four core sides, and $u_0=0$, $\Delta u=1$ for the
annulus sides. The windows belong to the topology. The interface map
$\xi_1(t)$ of (16), the restriction R and E therefore do not depend on the
curves. A new curve moves the side end vertices onto C, and on the O-grid
the inner corners follow as $c+\beta(v_{k+4}-c)$ with a fixed $\beta$. Only the
block maps are rebuilt: the annulus map is ruled between the periodic
Greville interpolants of wall and interface, the core interface sides take
the exact restriction (16) of the annulus interface curve, and the core
blocks are ruled between that side and the opposite straight side (a Coons
patch with straight radial sides).

Test (p = 3): D shapes $R=1.7+a\cos(\theta+\delta\sin\theta)$,
$Z=\kappa a\sin\theta$ with $\kappa=1.8$, $\delta=0.5$, wall a = 0.8,
interface a = 0.5, N = 256 samples. Neumann problem $u-\Delta u=f$ with exact
boundary flux on the spline wall:

| n | L2 error | area error |
|---|---|---|
| 4 | 1.8e-3 | 4.2e-3 |
| 8 | 3.6e-4 | 2.4e-4 |
| 16 | 8.1e-6 | 8.7e-6 |
| 32 | 2.8e-7 | 4.6e-7 |

The step 4 to 8 is slow (rate 2.3) because the error at n = 8 is dominated
by the annulus and the ring blocks next to the D tips. From n = 8 on the
rate is 5.5 and 4.9. After a new interface curve
($\kappa=1.6$, $\delta=0.3$) the space keeps its extraction, the maps follow
the new vertices exactly, and the wall area is unchanged to 1e-10.

## Free-boundary coupling

The vacuum outside the spline domain Omega with boundary Gamma enters through
the exterior Steklov-Poincare operator S of the axisymmetric Grad-Shafranov
operator $L\psi=-\nabla\cdot(R^{-1}\nabla\psi)$, discretised on quadratic
boundary panels with P2 traces (Costabel form
$S_h=W+(M/2-K)^TV^{-1}(M/2-K)$). The panel nodes $x_n$ are the images of the
cell-side ends and parameter mid points on Gamma, ordered counter-clockwise.
With the trace matrix $T_{nI}=\varphi_I(x_n)$ the free-boundary energy is

$$
\mathcal E(c)=\sum_{e,q}w_q\Bigl(\frac{|\nabla\psi|^2}{2R}-f\psi\Bigr)
+\tfrac12\,(Tc)^TS_h\,(Tc),
\tag{14}
$$

and its Hessian is $A+T^TS_hT$. The minimiser satisfies $L\psi=f$ in Omega
and $R^{-1}\partial_n\psi=-S\psi$ on Gamma, i.e. $\psi$ continues as a
decaying vacuum field. Test: $w=S_5(t)\,G(x,x_0)$ with the Green function G,
$t=(|x-x_0|-a)/(b-a)$ and the C5 smoothstep $S_5$, so $w=G$ outside a disc
inside Omega, and $f=Lw$. For p = 3 on the arc O-grid with 2 panels per cell
side, the L2 errors are 1.4e-2, 1.5e-3, 5.2e-5 for 4, 8, 16 cells per block
direction. Scaling $S_h$ by 1.1 stalls the error at 9e-2.

## Open

- Mixed regularity, $C^{p-1}$ inside blocks with $C^{p-2}$ only near C1
  interfaces, would keep the higher interior smoothness. The papers read
  assume uniform r.
- Vertex functions are scaled by $\sigma^{j_1+j_2}$ with the $\sigma$ of
  arXiv:1812.09088 (p. 13). The exponent is not legible in the extracted
  text. Any scaling spans the same space.
- THB on C1 spaces (needs the admissibility conditions of Bracco, Giannelli,
  Kapl, Vázquez, arXiv:2204.10000) and a built-in a posteriori estimator.
- Legs of a secondary X-point, C1 on X-point layouts.
- Preconditioning of THB matrices: Jacobi CG needs 705 iterations at 2686
  DOFs on six levels against 469 at 6229 uniform DOFs.
- Toroidal: sparse mode coupling for periodic splines, threading,
  stellarator-symmetric mode sets, field-period guide curves (nfp > 1 only
  for the cylindrical chart), section convergence on the figure-8 tube at
  higher M.
- Locality of the vertex functions at split vertices: the null-space basis
  (18) is orthonormal, not sparse.
- Corner parameters chosen from the curve (e.g. equal angle about the axis)
  instead of fixed quarters of the chord length.
- C1 up to a smooth wall without C0 lines needs a curved analysis-suitable G1
  reference (Kapl, Sangalli, Takacs, CAD 2018, not yet read).
- Layout generation from O-points, X-points and curves.
