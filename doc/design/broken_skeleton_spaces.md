---
title: Broken and skeleton space layouts
---

# Broken and skeleton space layouts

`fortfem_broken_skeleton_spaces` owns the discrete topology needed before a
broken finite-element or hybridized formulation can assemble. It deliberately
does not evaluate basis functions or choose a PDE. The same maps can therefore
be used by H¹, H(curl), H(div), L², DG, HDG, cut/XFEM, and IGA callers.

## Broken volume spaces

`initialize_broken_space_layout` creates one independent block of
`local_dof_count` degrees of freedom for every cell. Numbering is deterministic
and cell-major:

\[
  g(i,K)=(K-1)n_\mathrm{local}+i.
\]

The optional sign map records orientation of edge or face moments. H¹ and L²
layouts reject orientation flips; H(curl) and H(div) layouts accept the frozen
`+1/-1` signs needed by Piola and FEEC maps. The maps can be copied out with
`broken_space_layout_maps` for caller-owned local assembly.

## Skeleton trace spaces

`initialize_skeleton_space_layout` creates one shared trace block per facet and
records two adjacent cell IDs. A boundary side uses cell ID zero and sign zero;
an interior side uses a valid cell ID and sign `+1` or `-1`. The facet-major map

\[
  g(i,F)=(F-1)n_\mathrm{trace}+i
\]

is suitable for scalar, tangential, or normal traces. `skeleton_space_layout_maps`
exports the global facet map together with the side incidence and orientation,
which can then be passed to the dense reference or sparse CSC HDG assemblers.

## Scope and verification

The layout is a fixed-topology contract. Mesh connectivity, basis evaluation,
trace laws, metric transformations, static condensation, and sparse ownership
remain separate layers. The independent test checks deterministic numbering,
orientation preservation, boundary-side encoding, and rejection of invalid
dimensions and signs.
