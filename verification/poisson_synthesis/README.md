# FortFEM weak-form synthesis pilot (issue #62)

A self-contained first vertical slice for deriving FEM operators from
first-principles mathematical definitions, generating standard Fortran, and
verifying the algebraic/geometric identities that make the element kernel
correct.  This addresses GitHub issue #62 ("Pilot Fortran Synthesis for
weak-form derivation, verified element kernels, and FEM invariants").

The pilot has no external dependencies: it builds and runs with plain
`gfortran`/`fpm`.

## What it does

For the canonical 2D scalar Poisson problem on P1 (linear) triangles:

1. **Formal spec** — encodes the strong form `-Delta u = f`, the weak form
   `int grad(u).grad(v) = int f v`, the P1 barycentric basis, and the affine
   reference-to-physical map directly in the generator.
2. **Derivation** — derives the reference stiffness in closed form and verifies
   its patch-test properties (symmetry, zero row sums) before emitting anything.
3. **Generated kernel** — emits `src/generated/fortfem_poisson_p1_kernel.f90`,
   a dependency-free module with the element metric, stiffness, load, residual,
   and Jacobian, from the chain-rule physical gradients.
4. **Verified invariants** — the test checks partition of unity, `g^-1 g = I`,
   `det(g) = J^2`, symmetry of the bilinear form, residual/Jacobian consistency
   by finite differences, and the constant-field patch test.
5. **Source readback** — the committed `.f90` is re-read from disk and its
   kernel is compared against an independent adjugate-metric derivation, so the
   generated translation is validated against a second arithmetic path.
6. **Convergence** — a manufactured `sin(pi x) sin(pi y)` solution on the unit
   square converges at the P1 rates `O(h^2)` in L2 and `O(h)` in H1.

## Files

| File | Purpose |
|------|---------|
| `app/gen_poisson_p1_kernel.f90` | generator: spec, derivation, emission |
| `src/generated/fortfem_poisson_p1_kernel.f90` | generated standard-Fortran kernel |
| `test/test_poisson_p1_synthesis.f90` | invariant + readback + convergence tests |
| `fpm.toml` | dependency-free package manifest |

## Build, regenerate, test

```sh
cd verification/poisson_synthesis
fpm test          # runs the verification suite
fpm run --target gen_poisson_p1_kernel   # regenerate the kernel
```

The generator is deterministic; a regeneration produces a byte-identical
kernel.  `FORTSYNTH_OUTPUT_DIR` can override the output directory for
out-of-tree provenance checks.

## What the verification caught

During development the naive chain-rule implementation used
`grad_x phi = B^-T grad_xi phi`, which is correct only for axis-aligned right
triangles (where `B^-1 = B^-T`).  The independent adjugate-metric readback
check and the convergence suite both flagged the discrepancy on the 45-degree
corner triangles of the mesh.  The correct identity, with the map
`x = x0 + B^T xi`, is `grad_x phi = B^-1 grad_xi phi`.  This is exactly the
class of subtle geometry bug the issue's verification split is meant to catch.

## Trusted boundary and later scope

The formal verification layers named in issue #62 (FortSym/Lean identities,
Why3/SMT assembly contracts) remain aspirational and are not exercised by this
dependency-free pilot.  What is exercised here is a runnable, deterministic
synthesis pipeline: first-principles derivation -> generated Fortran ->
invariant/readback verification -> analytic convergence.  Extending the same
mechanism to nonlinear forms, DG fluxes, enrichment, IGA maps, exact-sequence
operators, and JVP/VJP generation (the issue's later-scope list) is
straightforward on top of this structure.
