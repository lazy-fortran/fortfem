## Spline and Fourier lines

`fortfem_bspline_lines` and `fortfem_fourier_lines` fill the shared
`tensor_line_tables_t` contract. Quadrature weights include span lengths and
`d` differentiates in the parameter coordinate. These are polynomial spline
and Fourier tables, with geometry and rational NURBS weights supplied by the
operator layer. Construction is setup work, outside the operator hot path.

`build_bspline_line(degree, line, status, ...)` accepts exactly one of
`nelem` (uniform spans on [0,1]) and `knots`. Open `knots` is a complete,
nondecreasing, clamped knot vector, including endpoint multiplicities p+1.
Interior multiplicities determine continuity and may range from 1 to p.
Periodic `knots` instead gives strictly increasing span breaks, including
both endpoints. `continuity` defaults to p-1 for generated knots and periodic
breaks and accepts 0 through p-1. If supplied for an open knot vector it
checks its interior multiplicities. Degree zero builds piecewise constants.
`nq` defaults to p+1 Gauss points per positive span. Invalid inputs return a
nonzero integer status. Nonuniform spans are supported in both cases.

The builder calls the existing `evaluate_bspline_basis` Cox-de Boor routine.
Periodic knots extend by translation, and local support maps to global dofs
modulo the number of periodic functions. Several local functions can map to
the same dof on a small periodic mesh, so consumers must scatter-add them.

`build_bspline_line_pair` requires p >= 1 and returns V0, V1 and a dense
coefficient incidence matrix D with shape (V1%ndof,V0%ndof). V1 removes the
first and last knot, lowers degree to p-1 and continuity by one (C0 becomes
discontinuous), and scales its functions as
`M_j = p B_{j,p-1}/(t_{j+p+1}-t_{j+1})`, with indices in the V0 knot vector.
Thus `(D u)_j = u_{j+1}-u_j`, wrapping at the periodic seam, and the
reconstructed derivative equals the V1 field at every point inside a span.
Open pairs have one fewer V1 dof. Periodic pairs retain equal dimensions and
a constant cokernel, as required by the topology of a circle. Tensor products
of the directional V0/V1 pairs provide the patch H1/H(curl)/H(div)/L2 spaces.

`build_fourier_line(M, line, status, nq, nfp)` uses one periodic cell on
[0,2*pi/nfp), ordered `1, cos(nfp*x), sin(nfp*x), ...`. Both nfp and nq are
optional and default to 1 and 2M+1, respectively. The equally weighted
trapezoidal rule requires nq >= 2M+1 and exactly integrates retained basis
products. Evaluation and derivatives reuse FortFEM's Fourier mode registry.
The mass diagonal is the period for the constant and half the period for
other modes, and stiffness divided by mass has eigenvalues (k*nfp)^2.
`build_fourier_line_pair` chooses V1 to be the derivative image, omitting
the constant. Each cos/sin block of D maps to `(-k*nfp*sin, k*nfp*cos)`.
For M=0, V1 has zero functions. This zero-mean V1 convention omits the harmonic
one-form of a full periodic de Rham complex. Consumers needing that harmonic
component should retain the constant using `build_fourier_line` for V1.

GLL derivative/histopolation pairs are not supplied here. The MF1 nodal line
builder provides V0, but a compatible degree-reduced edge basis remains a
separate extension.

`test_lab_spline_lines` checks partitions, polynomial moments, all basis
commuting identities, periodic wrapping and cardinal values, Fourier products
and eigenvalues, and 1D and 2D mass/stiffness against the existing sparse
B-spline assembly on an affine nonuniform patch, with a 1e-13 assembly tolerance.
