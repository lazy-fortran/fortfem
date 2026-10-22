# Laboratory problem contract

[H] All dimensions are nondimensional. The amplitudes below are model choices,
not reconstructions of an experimental equilibrium. No private inputs are used.

## Geometry (L1, L2, L3b, L4, L5)

[H] Let `(u,v)` range over the **closed unit disk**, and let `phi` be periodic
with period `2*pi`. Define

```
a(phi) = 1 + 0.1 cos(5 phi)
Rc(phi) = 3 + 0.1 cos(5 phi)
Zc(phi) = 0.1 sin(5 phi)
R = Rc(phi) + a(phi) u
Z = Zc(phi) + 1.3 v
X = (R cos(phi), R sin(phi), Z)
```

This defines the entire solid torus, including `(u,v)=(0,0)`, with five field
periods. It is an analytic shaped, phi-dependent stellarator-like test domain.
It is not a force-balanced stellarator. Its inverse is
`u=(R-Rc)/a`, `v=(Z-Zc)/1.3`. Define `s=u*u+v*v` without evaluating a poloidal
angle. All physical fields below are Cartesian or cylindrical functions.
The absolute Jacobian in `(u,v,phi)` is `1.3*a*R`, and `R>=1.8`.
The magnetic axis is the closed curve `(Rc(phi),phi,Zc(phi))`, not `R=0`.
A through-axis mesh contains that curve in its interior, with no duplicated
axis degrees of freedom. The affine tetra baseline approximates the boundary
by planar faces. Its norms are over the discrete domain `Omega_h`, so geometry
error must be qualified separately before claiming an error on the exact torus.

## L1: Poisson (criterion 1)

[H] Solve `-Delta U=f` with exact Dirichlet trace and toroidal periodicity, where
`U=sin(x)+cos(y)+sin(z)+(x+y+z)/10` and `f=sin(x)+cos(y)+sin(z)`.
The field is smooth through the magnetic axis and is independent of any flux
chart. This is an interior elliptic baseline relevant to FEM+BEM composition.
It does not by itself verify a BEM or DtN coupling.

## L2: H(curl) curl-curl plus mass (criterion 1)

[H] Use `A=(sin(y),sin(z),sin(x))`, `curl A=(-cos(z),-cos(x),-cos(y))`.
Then `div A=0` and `curl curl A + A = 2 A`. Prescribe its tangential trace
`n cross A` at the torus boundary, with periodicity in phi. Report relative L2
and relative H(curl) graph errors, including tangential trace lifting cost.
Nodal vector Poisson is not an admissible H(curl) baseline. The FortSym scalar
product returns the components of A as outputs 6:8, and the forcing is twice
those components. The adapter uses the existing first-order Nedelec curl-mass assembler,
edge-moment tangential trace lifting, and the covariant Piola map. Orders
above one are rejected until their face-moment lifting is supplied.

## L3a: two-dimensional Guenter/Sovinec closed-line test (criterion 3)

[V] The primary paper is S. Guenter, Q. Yu, J. Krueger and K. Lackner,
*Modelling of heat transport in magnetised plasmas using non-aligned
coordinates*, JCP **209** (2005), 354–370,
[doi:10.1016/j.jcp.2005.03.021](https://doi.org/10.1016/j.jcp.2005.03.021).
The body was checked in the author's published paper, pp. 357–361.
The cosine case used in its FEM comparison (p. 359 and Fig. 6) is

```
Omega = [-0.5,0.5]^2
T = psi = cos(pi*x) cos(pi*y)
B = ez cross grad(psi), b = B/|B|
K = I + (ratio-1) b b^T
-div(K grad T) = Q = 2*pi^2*T, T|boundary = 0
ratio in {1, 1e2, 1e4, 1e6}, chi_perp=1
```

At isolated magnetic nulls the direction is immaterial to the continuous
solution because `grad T=0`. The adapter uses `b=0` at a numerical null and
retains the isotropic tensor there. Quadrature does not place finite weight on
an isolated null. This convention must not regularize a neighbourhood of it.
The baseline uses FortFEM's existing `assemble_tensor_diffusion_matrix` on P1
triangles, with quadrature of the variable tensor, not the constant-tensor
example's right-hand side. It is expected to exhibit pollution.

[V] The paper also uses circular flux lines (p. 357, Fig. 2), centered at
`(0,0)` or `(-0.5,-0.5)`, with `T=1-r^3`, `Q=9*r` and exact Dirichlet data.
That case's central-temperature error is normalized by the exact temperature
drop to `(-0.5,0)`. It is a distinct case, not the cosine FEM comparison.

[H] For the cosine case, report signed
`chi_perp,num = chi_perp,eff - chi_perp = 1/T_h(0,0) - 1`, together with its
absolute value when comparing to Fig. 6's `|chi_eff-1|`. This uses the
amplitude-equivalent conductivity of this unit-amplitude source problem.
It is a diagnostic, not a tensor identification. Never replace L2 error with
this scalar or interpret a negative value as physical negative conductivity.
The 3D extension below uses a separate least-squares amplitude diagnostic.

## L3b: sheared three-dimensional field (criteria 2 and 3)

[H] Set `iota(s)=0.4+0.2*s` and `T=1-s`. In cylindrical orthonormal components,

```
t_R   = -0.5 sin(5 phi) (1+u) - a(phi) iota(s) v
t_phi = R
t_Z   = 0.5 cos(5 phi) + 1.3 iota(s) u
b = t/|t|
K = I + (ratio-1) b b^T
Q = -[T_RR + T_R/R + T_phiphi/R^2 + T_ZZ]
```

These definitions extend smoothly to the magnetic axis and satisfy
`t dot grad s=0`. They are the push-forward of
`d/dphi + iota(s)*(-v*d/du+u*d/dv)`. Thus `b dot grad T=0`, and Q is independent
of ratio. `gen_lab_problems` differentiates the inverse geometry to emit Q.
Prescribe the exact trace on the approximate mesh boundary and periodicity.
For a homogeneous exact-boundary solve, the amplitude diagnostic is
`alpha = integral(T_h*T)/integral(T*T)` and
`chi_perp,num = 1/alpha-1`. [H] This extends the 2D amplitude diagnostic and is
not the paper's 3D island transport diagnostic. For an approximate boundary,
report it explicitly as an amplitude proxy and include the trace error.

[H] An aligned comparison follows `theta(phi)=theta0+iota(s)*phi` within a
patch. At the periodic seam the radial-dependent shear generally fails to
close on the original nodes. A conforming interface interpolation or a
multi-patch constraint is required. A sinusoidally twisted mesh is **not** a
field-aligned mesh and is rejected by the current driver. Alignment is an
optional accelerator, never a condition on the physical problem.

## L4: nested map Newton model (criteria 2 and 4)

[H] Write `q=(s,theta,phi)` only for this comparison, with
`u=sqrt(s)*cos(theta)`, `v=sqrt(s)*sin(theta)`, and the same X as above.
Let `J_i = partial X/partial q_i`, `g_ij=J_i dot J_j` and
`G=sqrt(det g)`. For a physical map increment `w(q)` define

```
E(w) = 1/2 integral [C_ab g^ij partial_i w_a partial_j w_b
                     + lambda w_a w_a] G dq - integral F_a w_a G dq
C = [[2,1,0],[1,2,0],[0,0,1]], lambda=1
L(w)_a = -1/G partial_i(G g^ij C_ab partial_j w_b) + w_a
```

This is the exact Hessian of E. Its weak form has the curved metric, a
coupled SPD component block, a mass regularization and an imposed boundary,
which are representative costs of an elliptic regularized nested-map Newton
step. It omits the indefinite ideal-MHD force Hessian, gauge constraints,
pressure/current closures and nonlinear iteration, so it cannot predict
physical equilibrium convergence. C, lambda and the regularization are [H].
Use exact `w=(U,2U,3U)` and `F=C*(f,2f,3f)+w` from L1.
The Euclidean pull-back identity gives
`1/G partial_i(G g^ij partial_j (U o X)) = (Delta U) o X`.
This derives F independently of a chosen basis and retains the full metric,
not a diagonal-metric approximation. The scalar adapter diagonalises the
constant C exactly. Its eigenvalues are 3,1,1 and the squared exact modal
amplitudes are 4.5,0.5,9. It reuses the eigenvalue-1 solve for the two
proportional manufactured right-hand sides. Record three-component DOFs,
two unique assembled operators and their combined process time. This
specialised modal solve does not measure a variable-coefficient Newton block. In Cartesian meshes assemble the same
physical vector elliptic system. Polar-chart coefficients are singular at
s=0, whereas the underlying physical system is regular. Do not place a
quadrature point at s=0 or call this polar chart an axis solution.

[V] KIN6D context was read in `docs/equilibrium-nonnested-nn2.md`, including
its closed-line pollution and coupled-Newton discussion. The requested
`docs/island-equilibrium.md` and N1 handoff were absent at inspection time.
L4 does not import KIN6D implementation code or redefine its Newton equations.

## L5: axis regularity (criterion 2)

[H] Use L1's U, which has nonzero z gradient (`cos(z)+0.1`) and nonzero
curvature on the entire magnetic axis. Compare through-axis and polar-chart
families on the same physical problem. Report global L2 and H1 seminorms and
relative L2 error in `s<0.04` (a radius-0.2 tube in disk coordinates), with
physical volume weights. Also sample balls of physical radius 0.2 around axis
points if a local-ball comparison is needed. A tube is not a physical ball.
For every norm compute rates `log(e_coarse/e_fine)/log(h_coarse/h_fine)`.
Report rate loss as `rate_global-rate_axis` on each refinement pair and the
axis-rate difference between the two mesh families. Missing polar rows imply
no measured conclusion about chart-induced rate loss.

## L6: matched accuracy against VMEC-class solvers (criterion 4)

[V] Reuse [KIN6D](https://github.com/itpplasma/kin6d)'s
`benchmarks/nested/run_reference.py` and its README, without copying either
cases or adapter into FortFEM. The checked current KIN6D tree provides a
fixed-boundary Landreman–Paul QA case and integer-family oracle exports for
VMEC++, DESC and GVEC. The QA case has no exact equilibrium oracle. Its
integer-family comparator must never be applied to QA.

[H] Define matched accuracy by the **vector**
`(e_B,e_p,e_iota,e_flux)` with fixed common physical sample positions,
radial weights and normalizations. For integer-family cases, use the N1
exact oracle. For QA, use a jointly converged reference, retaining its
resolution and an independent refinement difference. Use Cartesian B and
physical pressure. Let `e_B=||B-Bref||_w/||Bref||_w`,
`e_p=||p-pref||_w/max(||pref||_w,Bref_rms^2/(2*mu0))`,
`e_iota=||iota-iotaref||_w/max(||iotaref||_w,1)` and
`e_flux=||psi-psiref||_w/|psi_edge-psi_axis|`.
Record flux sign/convention, boundary and profile hashes, physical sample
positions, convergence status, solve wall time, total wall time, actual
threads, peak process RSS, code version and compiler settings. A run matches
threshold epsilon only if **all four** components are <= epsilon and the
code reports convergence. Never substitute a force residual for e_B.
For zero-pressure QA, the magnetic-pressure denominator avoids division by
zero without silently dropping pressure. Align flux and iota conventions
before comparison. N1's present samples exclude axis and edge, so they do not
establish L5 regularity, volume norms or continuous error bounds.
Lab-side L6 equilibrium solves are a later KIN6D consumer milestone. No L1
or L4 solve is an L6 equilibrium timing.
