# Bíró tree--cotree curl--curl benchmark

This small manufactured case reproduces the structural part of the direct
tree--cotree gauge described by Bíró, Preis, and Richter, *On the use of the
magnetic vector potential in the nodal and edge finite element analysis of 3D
magnetostatic fields*, IEEE Transactions on Magnetics 32 (1996), 651--654
([DOI](https://doi.org/10.1109/20.497322)). It uses a cyclic graph-level
topology analogue, fixes the spanning-tree coefficients, forms a unit
topological curl--curl Gram block, solves the reduced cotree system directly,
and reconstructs the physical edge potential. The paper is a formulation
paper; its application-specific 3D geometry, metric, and source data are not
redistributed here.

The first output is `solution.png`: the graph and oriented edge-potential
arrows. `solution.csv` and `diagnostics.csv` are generated under
`output/example/biro_tree_cotree_benchmark`; no media are checked in. The
matrix and right-hand side are manufactured, so this is a reproducible
method/orientation oracle rather than a claim to ship the paper's source or
benchmark geometry.

Run with:

```text
fo run --example biro_tree_cotree_benchmark
```
