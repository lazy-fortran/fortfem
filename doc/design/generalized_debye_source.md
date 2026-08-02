---
title: Generalized-Debye-source coordinate contract
---

# Generalized-Debye-source coordinate contract

`fortfem_generalized_debye_source` is the topology- and kernel-neutral
coordinate layer for a BIEST-like vector surface formulation.  A tangential
coefficient vector is written as

\[
  j = G a + C b + H h,
\]

where `G` and `C` are caller-owned scalar-to-tangential lifts and `H` contains
one or more harmonic representatives of a closed surface.  The caller also
supplies a generalized-source map `S` and a period map `P`.  FortFEM returns

\[
  r_S = S j - s_\mathrm{target}, \qquad
  r_P = P j - p_\mathrm{target}.
\]

The source channels are deliberately not named as physical charge, current,
or magnetic flux.  Their units, differential operators, harmonic basis, and
normalization remain properties of the consuming BEM, FEM/BEM, DtN, or IGA
adapter.  This makes the same block usable on an RWG/BC pair, a high-order
surface space, or a spline patch graph.

The value routine and its JVP/VJP differentiate every supplied map, source
coefficient, harmonic coefficient, and target.  The VJP uses the real-part
complex pairing

\[
  \operatorname{Re}\sum_i \overline{\bar y_i}\,\dot y_i,
\]

and is checked by an independent manufactured dot-product oracle in
`test_generalized_debye_source_ad`.  This is a coordinate/lift layer only: a
second-kind Green-kernel equation, Beltrami parameter, resonance policy, and
tree--cotree or harmonic nullspace selection remain separate contracts.

The decomposition follows the generalized Debye-source strategy used by
BIEST-style integral formulations; see [Malhotra et al.](https://arxiv.org/abs/1902.01205)
and [O'Neil and Cerfon](https://arxiv.org/abs/1611.01420).  Those references
motivate the separation, but FortFEM does not implement a Taylor-state or
MRxMHD solver.

`assemble_generalized_debye_source_second_kind` is the matrix-free composition
boundary for a caller-owned second-kind surface block `K` (for example an
identity-plus-MFIE block):

\[
  r_K = K j - k_\mathrm{target}.
\]

It returns `r_K` together with the source and period residuals, and its JVP/VJP
propagate through `K`, all three lifts, both constraint maps, and every source
and target.  The block deliberately does not assemble a dense global matrix;
large BEM/IGA clients retain their own sparse, H-matrix, or FMM representation.
`test_generalized_debye_second_kind_ad` checks central reassembly and the full
real complex adjoint identity.
