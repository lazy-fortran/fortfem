# Linear-response interchange contract

`fortfem_linear_response_interchange` is the neutral data boundary for a
caller that supplies a linear ideal/resistive, vacuum, or wall response.  It
does not read GEQDSK/COCOS data, choose an equilibrium, or implement a GPEC,
MARS-F, GLISS, or STARWALL model.

For a fixed modal basis and topology, the public block composition is

\[
  L(\omega) = K - \omega^2 M + i\,\omega R + V + W,
\]

where `K`, `M`, and `R` are caller-owned equilibrium, inertia, and resistive
blocks, and `V` and `W` are caller-owned vacuum and conducting-wall blocks.
The frequency may be complex.  The module also carries integer `(m,n)` mode
pairs, labels, a response-channel matrix, provenance, and a positive
normalization scale.  It validates dimensions, finite values, duplicate modes,
and duplicate labels without imposing a physical normalization.

The JVP differentiates all five blocks and the frequency.  The VJP uses the
real part of the complex inner product,

\[
  \langle a,b\rangle = \operatorname{Re}\sum_i \overline{a_i}b_i,
\]

which is the convention used by the other complex FortFEM operators.  The
frequency adjoint is returned as a complex scalar under that convention.

The residual layer accepts a caller-owned state and source,

\[
  r(u) = L(\omega)u-b,
\]

and exposes its value, JVP, and VJP without owning a linear solver or a
stopping criterion.  This is the composition point for an external forced
response, eigen-residual linearization, or FEM/BEM/DtN/PML client.

`test_linear_response_interchange` is an independent manufactured oracle: it
checks the block signs, central differences, deep-copy semantics, duplicate
mode rejection, and the complex adjoint identity.  The contract is intended
to be composed with FEM, BEM, DtN, PML, Fourier, and tree--cotree blocks by an
external application.

`evaluate_linear_response_diagnostics` reports a normalized transpose
reciprocity defect and a conservative lower bound for the Hermitian part of a
response matrix.  A nonnegative bound certifies passivity for this Gershgorin
test; a negative bound is reported as inconclusive rather than being promoted
to a physical instability claim.  Exact reciprocity, passivity, and
normalization conventions remain properties of the external response model.

## Provenance

The block split follows the generic forced-response decomposition used by
linear ideal/resistive MHD and resistive-wall response codes.  FortFEM keeps
the implementation algebraic so that external adapters can record their own
equilibrium and code revisions.  Relevant parity targets are listed in
[`ROADMAP.md`](../../ROADMAP.md): GPEC, MARS-F, GLISS, and STARWALL.
