---
title: Fitted trace constraint
---

# Fitted trace constraint

`assemble_fitted_trace_constraint` assembles the multiplier block for a
duplicated fitted interface:

\[
 C=\begin{bmatrix}M(\lambda,T^+) & -M(\lambda,T^-)\end{bmatrix},\qquad
 M(\lambda,T)=\int_\Gamma \lambda\,T\,dS.
\]

The multiplier, plus, and minus trace bases may have independent column
counts. Only their quadrature rows and positive surface weights are shared.
The explicit sign is the oriented jump convention, so reversing the side
labels changes the constraint predictably instead of hiding a sign in a
global assembler. This is suitable for Lagrange multipliers, mortar traces,
fitted duplicated FEM spaces, DG/HDG skeletons, cut interfaces, and IGA
patches. Constitutive jump targets and physical boundary laws remain outside
the module.

Value, fixed-topology JVP, and real VJP actions cover all three trace bases
and the surface weights. The independent test checks a direct signed
cross-mass oracle, central differences, the dot-product identity, and
rejection of non-positive surface weights.
