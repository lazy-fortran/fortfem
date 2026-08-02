# Resistive-MHD branch-history contract

`resistive_mhd_branch_history_t` is a small, closure-neutral data object for a
continuation path.  A caller supplies strictly monotone continuation
parameters, state samples, residual samples, energy samples, and integer
branch labels.  The contract deliberately does not select a branch, run a
nonlinear solver, or prescribe a plasma model.

`evaluate_resistive_mhd_branch_diagnostics` counts distinct labels in one
history.  `compare_resistive_mhd_branch_histories` compares two histories on a
common fixed parameter grid, counts the union of labels, and reports maximum
state, residual, and energy differences.  The discrete multiplicity and
`hysteresis_detected` flag are diagnostics only; they are not differentiated.

For a continuous fixed-topology quantity,
`evaluate_resistive_mhd_branch_path_metric` integrates the caller-owned

\[
  f_i = E_i + \tfrac12(\|u_i\|^2+\|R_i\|^2)
\]

with the trapezoidal rule in the supplied continuation parameter.  Its JVP and
VJP implement the same rule analytically, including parameter, state,
residual, and energy directions.  This metric is a neutral path diagnostic,
not a physical energy functional; clients may replace it with their own
objective or ledger.

The independent test checks strict-monotonicity rejection, multiplicity and
hysteresis on two manufactured histories, central differences for the JVP,
and the full parameter/state/residual/energy dot-product identity for the VJP.
The primitive composes with the existing pseudo-arclength and deflated
residual contracts without duplicating their continuation or branch-selection
policies.
