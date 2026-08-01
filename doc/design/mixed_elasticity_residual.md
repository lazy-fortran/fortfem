---
title: Mixed elasticity residual
---

# Mixed elasticity residual

`assemble_mixed_elasticity_residual` is the neutral first-order algebraic
contract for mixed stress--displacement formulations. Given caller-owned maps

\[
 C : \sigma\mapsto\text{constitutive coordinates},\qquad
 E : u\mapsto\text{strain coordinates},\qquad
 D : \sigma\mapsto\text{equilibrium coordinates},
\]

it evaluates

\[
 r_c=C\sigma-Eu,\qquad r_e=D\sigma-f.
\]

`C` can be a compliance, a weak-symmetry block, or a strongly anisotropic
caller-owned constitutive matrix. `E` and `D` can be assembled from an
elasticity complex, TDNNS trace map, Hellinger--Reissner element, or an
external model. FortFEM does not choose the element family, symmetry law,
boundary condition, or material closure.

The JVP differentiates all six inputs. The VJP returns cotangents for both
residual blocks and all maps, including

\[
\bar\sigma=C^T\bar r_c+D^T\bar r_e,\qquad
\bar u=-E^T\bar r_c,\qquad \bar f=-\bar r_e.
\]

The focused test uses independent matrix-product expressions, a forward
finite difference of every map and state direction, and a complete real
adjoint identity. It also rejects incompatible dimensions. This makes the
block suitable for later compatible elasticity, tensor-pressure waves, and
acoustic/elastic interface compositions without claiming a full elasticity
solver.

The formulation is aligned with the mixed weak-symmetry construction of
[Arnold, Falk, and Winther](https://arxiv.org/abs/math/0701506), the
[JKU Linz TDNNS work](https://www.numa.uni-linz.ac.at/Talks/abstract/154/),
and the anisotropic mixed finite-element results of
[Pechstein and Schöberl](https://doi.org/10.1002/nme.3319).
