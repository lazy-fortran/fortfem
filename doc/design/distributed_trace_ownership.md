---
title: Distributed trace ownership and reduction
---

# Distributed trace ownership and reduction

`distributed_trace_layout_t` is the serial-first ledger for a physical trace
whose rows are split across several local partitions.  Each partition is a
`partition_layout_t`, so local-to-global row IDs, owner ranks, and explicit
owned/ghost masks are validated before any residual or mortar values are
reduced.  The ledger stores only fixed metadata and packed row offsets; it has
no communicator, mesh, geometry, or solver state.

The initialization contract requires that every global row appearing in the
ledger has one owner, and that all copies (including ghosts) name the same
owner rank.  A later MPI or task backend can exchange the packed local rows
without changing the physical FEM, BEM, DtN, DG, or IGA kernels.

## Value and derivative actions

For local rows (v_r) with component index (c),

\[
  V_{gc} = \sum_{r:g(r)=g} v_{rc},\qquad
  V^{\rm own}_{gc} = \sum_{r:g(r)=g,\;r\;\mathrm{owned}} v_{rc}.
\]

`assemble_distributed_trace_reduction` returns both arrays.  The first sums
duplicate ghost contributions in deterministic partition/row order.  The
second is the unique owner contribution and is useful for owner-only residual
rows or diagnostics.  The JVP applies the same fixed map to row increments.
The VJP takes cotangents for both outputs and adds the global cotangent to
every local copy, adding the owner cotangent only to an owned row.  Thus the
real transpose identity is explicit and independent of a communication layer.

Ownership, IDs, masks, row counts, and component counts are discrete fixed
metadata.  They are rejected when inconsistent rather than differentiated.
Physical coordinates, normals, surface metrics, trace basis values, and
constitutive/interface laws remain caller-owned and compose with this ledger.

`test_distributed_trace_ownership` builds two independent rank layouts,
checks an oracle that sums duplicate ghosts and owner rows, verifies a central
finite-difference JVP and the VJP dot-product identity, and rejects invalid
IDs, local masks, and cross-partition owner metadata.
