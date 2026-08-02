# Smooth deflated residual contract

`assemble_deflated_residual` is a geometry- and physics-independent wrapper
for clients that want to search for multiple branches of a nonlinear residual.
Given caller-owned residual samples `F(x)` and fixed reference states `z_j`, it
returns

\[
 F_d(x)=D(x)F(x),\qquad
 D(x)=1+s\sum_j(\lVert x-z_j\rVert^2+\delta^2)^{-p/2}.
\]

The positive `shift` \\(\delta\\) keeps this contract smooth at a reference
state. `scale` is non-negative and `power` is positive. Reference states are
fixed data in the derivative contract: changing or adding roots is a discrete
client policy and must not be differentiated.

The JVP accepts the caller's base-residual directional derivative and adds the
deflation-factor derivative. The VJP returns cotangents for the base residual
and state, which makes the wrapper usable in Newton--Krylov, adjoint, and
DESC-like objective paths without selecting a nonlinear solver. Combine it
with `assemble_pseudo_arclength_residual` when a continuation row is needed.

The focused test checks the multiplicative value formula, an independent
central-difference JVP oracle, the real adjoint identity, and invalid zero-shift
rejection. No plasma closure, branch-selection rule, or application file
format is implemented here.
