# Multi-region Beltrami residual contract

`fortfem_beltrami_residual` is a physics-neutral residual product for clients
that represent a relaxed magnetic region by a supplied curl eigenvalue
relation.  For region `r`, sample `q`, and Cartesian component `c`, it returns

\[
 R_{rqc}=({\rm curl}\,B)_{rqc}-\lambda_r B_{rqc}.
\]

The module does not evaluate a curl, choose a gauge, construct a vector
potential, or solve for a field.  A compatible finite-element, IGA, Fourier,
or external application layer owns those operations and supplies the sampled
`curl_field`, `magnetic_field`, and region scalar `lambda` arrays.

## Constraint rows

Optional rows are represented by value/target pairs.  The divergence pair is a
rank-2 array (arbitrary row and quadrature layout); flux and helicity are
rank-1 arrays.  Their residuals are simply

\[
 r_{\nabla\cdot B}=d-d_\star,\qquad
 r_\Phi=\Phi-\Phi_\star,\qquad
 r_H=H-H_\star.
\]

A zero-sized pair disables that family, so the same contract covers a
single-region solve, a multi-region solve, or a client that supplies only the
curl-eigen block.  This deliberately keeps flux and helicity definitions
external: region topology, loop representatives, normalization, and gauge
choices belong to the consuming application.

## Differentiation and validation

The JVP includes the product rule,

\[
 \dot R=\operatorname{curl}\dot B-dot\lambda B-lambda\dot B,
\]

and the VJP is the exact real transpose action.  Constraint row derivatives
are value-minus-target.  Every value and direction is checked for compatible
rank/shape and finite entries before any output is accepted; invalid calls
return `FORTSPARSE_INVALID_MATRIX` and zero outputs.

## Provenance and scope

The curl-eigen relation is the neutral algebraic block used by multi-region
relaxed-MHD formulations (for example, SPEC's Beltrami regions).  See the
[SPEC documentation](https://princetonuniversity.github.io/SPEC/) and the
MRxMHD formulation by Hudson *et al.*
([doi:10.1063/1.4977886](https://doi.org/10.1063/1.4977886)) for the physical
context.  FortFEM only provides the differentiable array contract; it does
not implement SPEC/SPECTRE physics, equilibrium profiles, readers, gauges, or
  solvers.

## Two-path manufactured parity

`compare_beltrami_two_region_residual` composes the residual above with an
independent algebraic oracle. A compatible H(curl) client supplies
`curl_hcurl`; the report's oracle path evaluates `curl_oracle-lambda*B`
directly, so a shared residual implementation cannot make the comparison pass
by construction. The report records absolute and relative path errors and the
divergence, flux, and helicity residual norms.

`validate_beltrami_resonance` takes a caller-supplied table of forbidden
gauge/eigenvalue parameters and returns `FORTSPARSE_SINGULAR` when a region
parameter falls within the declared tolerance. It does not infer a spectrum
or choose a gauge, leaving that policy to compatible H(curl), BIEST-like,
Fourier, and IGA clients.

## Slab and toroidal-shell ledger

`compare_beltrami_shell_residual` is a small geometry-labelled acceptance
fixture for two-region slab and toroidal-shell clients. The label identifies
the caller's fixed physical sample set; FortFEM does not construct coordinates
or a shell mesh. It independently compares weighted
`curl_hcurl-lambda*B` and `curl_oracle-lambda*B` norms and closes supplied
flux, helicity, and per-region quadratic field-energy targets. The energy
quadrature is

\[
 E_r=\tfrac12\sum_q w_{rq}\,B_{rq}\cdot B_{rq}.
\]

The report marks `ledger_closed` only when all three target families close
within the declared tolerance. This is a neutral parity and conservation
oracle for BIEST-like, compatible H(curl), Fourier, and IGA implementations;
it does not choose an equilibrium, topology, or gauge.
