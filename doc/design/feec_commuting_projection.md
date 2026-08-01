---
title: FEEC commuting-projection diagnostics
---

# FEEC commuting-projection diagnostics

`assemble_feec_commuting_projection` audits the commuting diagram used by a
discrete de Rham construction. It is deliberately a neutral coefficient-space
contract: it does not assume a mesh, quadrature rule, constitutive tensor, or
application-specific field names. The same routine can therefore be used for
fitted finite elements, cut/XFEM spaces, IGA and multipatch maps, and periodic
reductions.

For discrete maps \(G_d,C_d,D_d\), continuous maps \(G_c,C_c,D_c\), and
projections \(P_0,P_1,P_2,P_3\), the returned defects are

\[
  E_0 = G_dP_0-P_1G_c,\qquad
  E_1 = C_dP_1-P_2C_c,\qquad
  E_2 = D_dP_2-P_3D_c.
\]

The maps are supplied as rectangular dense arrays so that dimensions make the
domain and codomain of every arrow explicit. Incompatible dimensions, empty
maps, and non-finite values are rejected before any product is formed. The
routine does not silently project or pad a map.

## Differentiation contract

`assemble_feec_commuting_projection_jvp` applies the complete product rule to
all ten maps. `assemble_feec_commuting_projection_vjp` is its real adjoint;
the focused test checks

\[
 \langle \bar E,\dot E\rangle
 = \sum_X\langle\bar X,\dot X\rangle
\]

with an independently chosen perturbation. A central-difference check of all
three defects is included as a second behavioral oracle. As with the other
topology contracts, the discrete map graph and its dimensions are fixed during
differentiation; a changed cut topology, refinement, or orientation is a new
discrete problem and must be rebuilt.

## Enriched and cut FEEC use

The diagnostic is not a claim that an arbitrary scalar enrichment is a
compatible vector enrichment. Clients should run it for each enriched scalar,
\(H(\mathrm{curl})\), \(H(\mathrm{div})\), and \(L^2\) space. A nonzero defect
is useful information: it either identifies an implementation error or records
the exact-sequence identity intentionally relaxed by a physical jump or
singular interface law. Trace regularity, jump terms, and distributional
sources remain separate contracts.

The implementation is a short matrix composition whose product-rule algebra is
equivalent to the expression FortSym emits for a generated map. Keeping the
audit independent of a particular generated kernel gives an oracle for
FortSym-generated, hand-written, and external-client maps alike.
