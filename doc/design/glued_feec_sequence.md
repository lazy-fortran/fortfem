---
title: Signed glued FEEC sequence contract
---

# Signed glued FEEC sequence contract

`assemble_glued_feec_sequence` is the dense reference composition between
cell-local compatible maps and a global numbering. The four integer maps use
signed IDs: the absolute value selects the global degree of freedom and the
sign carries the caller's orientation. Contributions from shared IDs are
summed before the two composition diagnostics are formed:

\[
G = \operatorname{glue}(G_e),\quad
C = \operatorname{glue}(C_e),\quad
D = \operatorname{glue}(D_e),
\qquad CG,\ DC.
\]

This is deliberately a dense, deterministic reference path. It is suitable
for small regression fixtures and for validating a future sparse assembler;
it does not own mesh connectivity, metric Piola maps, trace penalties, or
material blocks. Distinct IDs produce broken DG/HDG spaces, while shared IDs
and signed orientations produce conforming, cut, or multipatch spaces.

The JVP differentiates the local matrices at fixed integer maps. The VJP
first accumulates cotangents through `C*G` and `D*C`, then scatters them back
with the same signs. Topology/map changes are discrete events and are not
differentiated.

## API

```fortran
call assemble_glued_feec_sequence( &
    local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
    hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, &
    gradient, curl, divergence, curl_gradient, divergence_curl, status)
call assemble_glued_feec_sequence_jvp(...)
call assemble_glued_feec_sequence_vjp(...)
```

All local maps must have positive-size dimensions and all signed IDs must be
nonzero and within their corresponding global count. The focused test uses
two cells with a shared vertex and an orientation-reversed edge, compares all
outputs against an independent loop oracle, checks a central-difference JVP,
and verifies the real VJP dot-product identity.
