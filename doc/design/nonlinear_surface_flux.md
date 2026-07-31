---
title: Nonlinear material-surface flux
---

# Nonlinear material-surface flux

`assemble_nonlinear_surface_flux` is the application boundary for wall,
target, sheath, recycling, radiation, and conjugate-transfer laws. FortFEM
does not encode any of those physical models. An application supplies a local
trace callback

\[
  F_q = F(u_q,n_q,\mathrm{tag}_q),
\]

while the library performs the common oriented quadrature and finite-element
pairing,

\[
 R_{ia} \mathrel{+}= T_{qia}\,w_q F_{qa}, \qquad
 L_{\tau a} \mathrel{+}= w_q F_{qa}.
\]

`L` is an inspectable integrated ledger for each positive integer surface tag.
The same tag and normal orientation are used in the residual, so a caller can
combine the ledger with an open-line FCI balance without introducing a ghost
plane or a smeared volume source.

The JVP callback supplies the local derivative of `F` with respect to the
trace state and normal. The assembly differentiates trace bases and
quadrature weights as well, using fixed tags. The VJP callback returns local
state and normal cotangents; the assembly adds the residual and ledger bars and
returns the complete real dot-product adjoint. Tag changes are discrete
topology events and are rejected by the shape/status checks rather than
silently differentiated.

The focused test uses two nonlinear tagged laws, an independent value oracle,
central differences for the complete JVP, and a full residual-plus-ledger
dot-product identity for the VJP. Invalid weights, missing tags, incompatible
trace shapes, non-finite data, and callback failures are rejected.

This contract is the surface-law complement to the PARALLAX-aligned FCI
geometry and terminal-flux primitives in
[`fci_parallel_operator`](fci_parallel_operator.html). It deliberately does
not implement a sheath model, particle pusher, reaction network, or material
database.
