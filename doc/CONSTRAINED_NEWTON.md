# Constrained Newton and parameter continuation

Two solver-orchestration modules for nonlinear clients that have an energy
with exact derivatives. They own no physics and no discretisation.

| Module | Content | Test |
|---|---|---|
| `fortfem_constrained_newton` (`src/solvers/constrained_newton.f90`) | Newton for a stationary point of W(z) under linear equality constraints and fixed slots, by LDLᵀ of the KKT matrix with its inertia; three step acceptances | `test/test_constrained_newton.f90` |
| `fortfem_parameter_continuation` (`src/solvers/parameter_continuation.f90`) | natural-parameter continuation with step control for a family P(s), s from 0 to 1 | `test/test_parameter_continuation.f90` |

## KKT step and inertia

A step solves, on the free slots,

    [ H + σ d I   Cᵀ ] [ δz ]   [ −g      ]
    [ C           0  ] [ λ  ] = [ d − C z ]

with H the Hessian (the caller gives its upper triangle), C the rows of the
constraints C z = d, d the smallest Hessian diagonal and σ a relative shift.
The matrix is factorised by the symmetric LDLᵀ of `fortfem_sparse_direct`
(FortSparse, Bunch–Kaufman pivoting), so zero diagonal blocks and indefinite
Hessians are admitted and the inertia is a by-product. With m independent
constraints the KKT matrix has m more negative eigenvalues than the Hessian
restricted to the null space of C (Sylvester's law applied to the block
factorisation; e.g. Nocedal and Wright, *Numerical Optimization*, 2nd ed.,
Theorem 16.3, cited from memory). The report gives that difference at the converged state, from
a factorisation without shift: the number of unstable directions of the
stationary point on the constraints.

Constraints for a continuous symmetry. When W is invariant under a group, the
tangents of the group orbit are null vectors of H at every stationary point.
A constraint row that is not orthogonal to such a tangent (a moment of the
update with the generator) makes the KKT matrix regular without a penalty
weight; its multiplier vanishes at the solution because the symmetry costs
no energy. A row that vanishes on the free slots is inactive (its multiplier
is kept at zero and it does not count in the inertia), so a client can keep
a fixed pattern when a generator degenerates.

Constraints that follow the state. With the optional `update` procedure the
rows are rebuilt at every iterate. With d = C(z_k) z_k every step lies in the
null space of C(z_k): the step is orthogonal to the generators at the current
state, and no reference state is carried along. The limit satisfies
g + C(z)ᵀ λ = 0 with the rows at the limit. The pattern of the rows must not
change.

Convergence: the gradient of the Lagrangian g + Cᵀλ, with the multipliers of
the step, scaled by the square roots of the Hessian diagonal, below the
tolerance, and the constraints satisfied.

## Step acceptance

| Option | Rule | Use |
|---|---|---|
| `NEWTON_ACCEPT_ENERGY` | the full step when W falls by 1e-4 of the fall of the quadratic model; otherwise σ × 10 (σ / 100 after an accepted step) | minima |
| `NEWTON_ACCEPT_ENERGY_RELAXED` | as above, with relaxed steps (below) | minima in a curved valley of a stiff term |
| `NEWTON_ACCEPT_GRADIENT` | the step a δz, a = 1, 1/2, …, when the scaled gradient of the Lagrangian falls below (1 − a/4) of its value | the stationary point next to the start, any signature |

Relaxed steps (the watchdog technique: Chamberlain, Powell, Lemaréchal and
Pedersen, *The watchdog technique for forcing convergence in algorithms for
constrained optimization*, Math. Programming Study 16 (1982); cited from
memory, to be confirmed). A Newton step along a curved valley of a stiff
term κ c(z)²/2 leaves the valley by the curvature and pays κ times the
fourth power of its length, so the energy test alone forces steps that
shrink with κ, although the next step returns to the valley. A step that
fails the energy test while the matrix is positive on the constraints (known
from the inertia, at no cost) is therefore taken on trial and kept when the
step after it brings W below the value required of the first. Otherwise the
state returns to the start of the pair, σ is raised, and no further trial is
made for a number of accepted steps that doubles with every failed pair, so
the wasted work stays bounded. Rejections on an indefinite matrix never start
a trial.

Measured (`test_constrained_newton`, W = κ (|u|² − 1)²/2 − x in the plane,
start at 108° on the circle, tolerance 1e-9):

| κ | steps, energy test | steps, relaxed | accepted pairs |
|---|---|---|---|
| 1e2 | 20 | 17 | 3 |
| 1e4 | 104 | 37 | 12 |
| 1e6 | 568 | 109 | 36 |

Both grow with κ; the relaxed test needs a third to a fifth of the steps.
Two alternatives were tried on the same case and are not in the code:

- the natural monotonicity test (Deuflhard, *Newton Methods for Nonlinear
  Problems*, Springer 2004, cited from memory: the simplified Newton correction
  from the factors of the step must contract) as a fallback of the energy
  test accepted no step that the energy test refused, at any κ;
- a second-order correction (the simplified Newton correction added to the
  trial point, then the energy test again) likewise.

Both use the factors of the old point. By an estimate on this case (not
measured separately) they are valid only while the stiff direction has
turned by less than about κ^(−1/2) along the step, which is not longer than
what the energy test admits. Neither is a remedy for a stiff
penalty term; continuation in the stiff parameter is.

## Parameter continuation

`continuation_run(stage, x, options, report, status)`: the caller's `stage`
solves P(s) from a given state and returns a status; the driver owns the
sequence of s. P(0) is solved from the start (with `try_target_first` P(1) is
tried before and P(0) only when that fails). From the last converged state at
s_ok the stage s_ok + ds is solved, ds = `max_step` at first; a converged
stage doubles ds up to `max_step`, a failed one halves it and is retried
from the same state; a failure with ds below `min_step` ends the run with the
last converged state. The end state solves P(1) whatever the path. No stage
is set by hand.

`test_parameter_continuation`: atan(K(s) (x − 1 − s)) = 0 with K = kmax^s,
whose undamped Newton basin is |x − 1 − s| < 1.39/K. For kmax = 10^1.5 the
direct solve at s = 1 fails and the driver arrives in 28 stages, 15 of them
failed (the doubling after every success is paid by a failure when the basin
shrinks monotonically); for kmax = 1e4 the basin falls below the smallest
step and the run stalls at the last converged state.
