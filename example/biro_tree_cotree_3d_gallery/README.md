# Bíró tree--cotree 3-D gallery

This executable is a small, license-safe manufactured 3-D edge-element
analogue of the direct tree--cotree gauge described by Bíró, Preis, and
Richter, *On the use of the magnetic vector potential in the nodal and edge
finite element analysis of 3D magnetostatic fields*, IEEE Transactions on
Magnetics 32 (1996), 651--654
([DOI](https://doi.org/10.1109/20.497322)).

The unit cube supplies the graph topology and six oriented face boundaries.
Edge unknowns are line integrals of a smooth manufactured vector potential;
the topological curl--curl Gram matrix is singular because discrete gradients
are invisible to curl.  FortFEM's spanning-tree restriction fixes seven tree
edge degrees of freedom, solves the five-dimensional cotree block directly,
and prolongs the result back to all twelve oriented edges.

The paper is a formulation paper.  Its application-specific 3-D geometry,
metric, material data, and source data are not redistributed here, so this
fixture should not be read as a claim to reproduce a paper figure.  It is a
reproducible method and orientation oracle that exercises the same direct
gauge mechanism in three dimensions.  The first plot is the solved physical
edge potential: blue tree edges, red cotree edges, black edge-potential
vectors, and edge-DOF magnitude markers.  The oriented face circulations are
reported separately in `diagnostics_1d.png` and `diagnostics.csv`.

Run it with:

```text
fo run --example biro_tree_cotree_3d_gallery
```

Generated files are written below
`output/example/biro_tree_cotree_3d_gallery`; media are intentionally not
checked in.
