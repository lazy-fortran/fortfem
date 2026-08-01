---
title: Elasticity weak-symmetry constraint
---

# Elasticity weak-symmetry constraint

`assemble_elasticity_symmetry_constraint` contracts a caller-owned fixed
topology map with a stress vector:

\[
 r_s=W\sigma-g.
\]

The rows of `W` may extract skew stress components for a weak-symmetry
multiplier, normal-normal traces for a TDNNS formulation, or another
application-owned symmetry constraint. The map is deliberately not tied to a
mesh, an element family, or a constitutive law.

The JVP is

\[
 Dr_s=\dot W\sigma+W\dot\sigma-\dot g,
\]

and the real VJP returns

\[
 \bar W=\bar r_s\sigma^T,\qquad
 \bar\sigma=W^T\bar r_s,\qquad
 \bar g=-\bar r_s.
\]

The focused test has an independent matrix oracle, forward-difference check,
complete adjoint identity, and dimension rejection. This block composes with
the mixed elasticity constitutive/equilibrium residual and with compatible
IGA or finite-element trace maps.

The contract follows the weak-symmetry perspective of
[Arnold, Falk, and Winther](https://arxiv.org/abs/math/0701506) and is
compatible with the TDNNS direction documented by the
[JKU Linz group](https://www.numa.uni-linz.ac.at/Talks/abstract/154/).
