# Enriched FEEC sequence contract

`assemble_enriched_feec_sequence` is the neutral matrix layer between shifted
XFEM/XIGA basis constructors and a client-owned FEEC assembly.  It accepts
base incidence maps

\[
G:V_0\to V_1,\qquad C:V_1\to V_2,\qquad D:V_2\to V_3,
\]

and four coefficient maps from enriched spaces into their base spaces:
`scalar_map` (`V0_e -> V0`), `hcurl_map` (`V1_e -> V1`), `hdiv_map`
(`V2_e -> V2`), and `l2_map` (`V3_e -> V3`).  The composed maps are

\[
G_e=V_1^TGS_0,\qquad C_e=V_2^T C V_1,\qquad D_e=V_3^T D V_2.
\]

The contract reports `C_e G_e` and `D_e C_e`; it does not silently project or
repair a non-commuting enrichment.  Thus shifted scalar and componentwise
shifted vector maps can be used with fitted, cut, Piola, or IGA constructions
without selecting geometry, global numbering, quadrature, or a physical PDE.

The JVP differentiates all seven caller-owned factors by the product rule.  The
VJP includes the two defect products and satisfies the real Frobenius
dot-product identity.  All dimensions are explicit and incompatible maps,
outputs, increments, or cotangents are rejected with `FORTSPARSE_INVALID_MATRIX`.
The implementation is a dense reference oracle; CSC or FortSym-generated
actions may implement the same algebra in a storage-specific layer.

The focused `test_enriched_feec_sequence` fixture compares value and JVP
outputs with independent dense `matmul` expressions, checks the VJP identity,
and exercises dimension rejection.  No physical PDE or global ownership is
hidden in this contract.
