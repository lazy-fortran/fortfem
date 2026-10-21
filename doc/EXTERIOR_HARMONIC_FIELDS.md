# Decaying exterior harmonic fields

`fortfem_exterior_harmonic_fields` constructs finite-energy representatives
on the complement of a connected genus-g handlebody. The caller supplies g
oriented, closed polygonal wires inside the excluded region. Their homology
classes must be independent; the wires must stay away from the exterior.
The constructor supports arbitrary curved boundary maps and P1/P2 traces.
Geometry, interior placement and linking independence are caller contracts.

With unit mu0 times current, each field and vector potential are

    H(x) = (1/4pi) integral_wire dl cross (x-y)/|x-y|^3,
    A(x) = (1/4pi) integral_wire dl/|x-y|.

Straight segments are integrated analytically. The generated leaves include
value/JVP/VJP. FortSym checks divergence, curl A = H, the total-derivative
identity whose endpoint contributions cancel on a closed wire, and the two
segment primitives. These identities hold away from the wires. Closed wires
cancel the monopole contribution: H = O(|x|^-3), A = O(|x|^-2).

The period on a boundary cycle is its oriented linking number with the wire.
The exterior has g independent periods, although the boundary has 2g cycles.
`prescribe_exterior_periods` solves on g selected independent cycles and
checks every supplied period. Status 3 rejects incompatible data. In
particular, on an unknotted torus an unlinked longitude has zero period for a
decaying exterior field; an arbitrary nonzero toroidal period cannot be
created by adding grad(phi), which is singular on the exterior axis.

`assemble_exterior_harmonics` checks closed oriented connected topology and
that the number of wires equals the genus. It returns the signed normal load
N and the untruncated exterior Gram matrix G:

    N_ij = integral_Gamma trace_i H_j.n,
    C = -N,
    G_ij = -integral_Gamma n.(A_i cross H_j).

The last identity follows from div(A_i cross H_j) = H_i.H_j in the exterior;
the sphere-at-infinity term vanishes by decay. Raw Gram asymmetry is returned
before averaging. Independently refine boundary quadrature and wire geometry;
this routine does not supply a certified error bound.

With the positive scalar exterior Steklov operator S,
`tangent_exterior_harmonics` solves S w = N and returns

    H_t,j = H_j + grad(w_j),
    G_t = G - N^T S^-1 N.

These representatives have zero normal trace in the discrete weak equation
and are orthogonal to scalar exterior gradients. Scalar corrections preserve
all periods. The physical normal-trace error also includes BEM discretization
and quadrature error; a small linear-system residual alone does not bound it.
Use the returned boundary trace w with the scalar BEM field reconstruction.
The API does not currently reconstruct the corrected field at exterior points.

`couple_harmonic_sector` consumes S, C, G in the same coordinates. For total
cut-flux constraints, include the scalar field's cut flux in the prescribed
right-hand side; a period constraint needs no such correction.

The independent tests check circular-loop axial fields, poloidal linking,
unlinked longitude, two handle periods, dipole decay, segment derivatives,
and boundary energy against an exterior volume integral using elliptic fields.
The volume oracle reports both quadrature refinement and a larger truncation
box. These are numerical accuracy checks, not interval certificates.

## Native validation checkpoint

At source 024ebeaad97be1c3a13144b96490117c61d51707, GNU Fortran 12.2
release and checked-debug runs pass all 54 assertions in the exterior field,
exterior energy, harmonic-energy and curved P2/P1 programs. The independent
checks include circular-loop axis values, linking periods, two handles,
dipole decay and boundary energy against elliptic-field volume quadrature.
The retained FortSym generator proves the off-wire identities and reproduces
the emitted leaf byte for byte.

The separate tangent-field reconstruction gives a relative normal-field L2
error of 0.00911163, below the unchanged 0.03 gate. Its corrected Gram value
is 2.359031529, versus 2.483421317 before the scalar correction, as required
by energy minimization. These are numerical observations for the retained
curved-torus fixture. They do not certify arbitrary geometry or establish
assembled geometry derivatives. Strong-scaling measurements remain in
progress; one eight-thread timing is insufficient to establish a speedup.
