---
title: Symplectic map defect contract
---

# Symplectic map defect contract

`assemble_symplectic_map_defect` is a neutral structure diagnostic for any
linear one-step map. Given a state map \(S\) and caller-owned form
\(\Omega\), it returns

\[
    E = S^T\Omega S-\Omega.
\]

An ideal Hamiltonian or symplectic update has \(E=0\) to the requested
floating-point tolerance. The routine does not assume canonical ordering,
choose a PDE, or classify a dissipative block as symplectic. This makes the
same check usable for mixed acoustics, vibration, elasticity, and compatible
electromagnetic states.

The JVP includes map and form directions:

\[
\dot E=\dot S^T\Omega S+S^T\dot\Omega S+S^T\Omega\dot S-\dot\Omega.
\]

The VJP uses the real Frobenius product and returns both map and form
cotangents. Fixed dimensions and finite values are required; changing the
state topology is outside this differentiable contract.

## API

```fortran
call assemble_symplectic_map_defect(map, symplectic_form, defect, status)
call assemble_symplectic_map_defect_jvp( &
    map, symplectic_form, map_dot, symplectic_form_dot, defect_dot, status)
call assemble_symplectic_map_defect_vjp( &
    map, symplectic_form, defect_bar, map_bar, symplectic_form_bar, status)
```

`test_symplectic_map_defect` compares the primal and tangent products against
independent matrix expressions, a central difference, and the real adjoint
identity. It intentionally uses a non-symplectic map so the defect is a
measured diagnostic rather than a tautological zero test.
