---
title: Partition layout and serial reduction
---

# Partition layout and serial reduction

`partition_layout_t` is the communicator-free data boundary for future
distributed assembly. It records a fixed local-to-global ID list, the owner
rank of every local entry, the local rank, and an explicit owned/ghost mask.
Local IDs are unique, owned entries must name the local rank, and ghost
entries must name another nonnegative rank. No MPI object is stored, so the
same local kernels run on a single node and inside a later halo backend.

`assemble_partitioned_sum` accumulates local values into a global reduction
buffer in local-index order. The result is deterministic and is the serial
no-op communicator implementation; a distributed client can replace the
reduction boundary with an all-reduce without changing element kernels.
The matching JVP and transpose VJP are linear map actions and do not
differentiate ownership or topology decisions.

The independent test checks owner/ghost metadata and exported maps, global-ID
accumulation, finite-difference JVP agreement, the real dot-product identity,
and duplicate-ID rejection.
