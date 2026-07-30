# Maxwell PML moving-mesh adjoint

This example benchmarks analytical forward and reverse products for
arbitrary-order-compatible tetrahedral Nedelec PML assembly. The differentiated
inputs include every mesh coordinate, each element's complex Cartesian
stretch, and the wave number.

The JVP propagates one simultaneous mesh/PML perturbation. The VJP consumes a
merged CSC matrix cotangent and returns the full design gradient in one reverse
element sweep, including accumulation at shared vertices. A complete
two-sided reassembly is retained as an independent directional oracle.

Generated gallery artifacts:

- `derivative_timing_1d.png`: measured JVP, VJP, and directional
  finite-difference timings, plus the coordinate-wise finite-difference
  full-gradient cost estimated from measured primal assembly;
- `adjoint_accuracy_1d.png`: relative reverse/central-difference error;
- `mesh_gradient_2d.png`: projected mesh-vertex descent directions;
- `benchmark.txt`: machine-readable sizes, timings, and errors.

No rendered image is committed.

## Method choice

Analytical reverse mode is the production choice for many geometry variables:
its element-sweep count is independent of the number of vertices. FortSym
generates the complex PML coefficient products, while FortNum's generated
determinant/inverse products and analytical Piola products handle geometry.
This path works with GCC. Enzyme remains an optional comparison oracle rather
than a runtime dependency.

## Provenance

- P. Monk, *Finite Element Methods for Maxwell's Equations*, Oxford
  University Press (2003), covariant Piola maps and Nedelec elements.
- J.-M. Jin, *The Finite Element Method in Electromagnetics*, 3rd ed.,
  Wiley (2014), Cartesian perfectly matched layers.
- J. Nocedal and S. Wright, *Numerical Optimization*, 2nd ed., Springer
  (2006), adjoint reduced gradients.
