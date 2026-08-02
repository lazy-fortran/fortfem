---
title: Signed glued FEEC CSC contract
---

# Signed glued FEEC CSC contract

`assemble_glued_feec_sequence_csc` is the sparse companion to the dense
`assemble_glued_feec_sequence` reference path. It accumulates the three local
de Rham maps into rectangular FortSparse CSC matrices using signed
local-to-global IDs. Duplicate entries from shared cells are compressed by
`csc_from_triplet`; the dense path remains the independent behavioral oracle.

The JVP reassembles only the local matrix directions at fixed topology. The
VJP reads canonical CSC cotangents and scatters them to local matrices using
the same row/column orientation signs. No sparse pattern or graph rebuild is
differentiated.

## API

```fortran
call assemble_glued_feec_sequence_csc( &
    local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
    hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, &
    gradient, curl, divergence, status)
call assemble_glued_feec_sequence_csc_jvp(...)
call assemble_glued_feec_sequence_csc_vjp(...)
```

The focused CSC fixture compares every matrix action against an independently
loop-assembled dense matrix, checks the fixed-topology tangent, verifies the
cotangent scatter, and rejects zero IDs. This is a small deterministic
reference for a future distributed sparse owner/ghost assembler; it does not
select a solver or PDE.
