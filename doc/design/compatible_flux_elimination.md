---
title: Compatible flux elimination
---

# Compatible flux elimination

`assemble_compatible_flux_elimination` is the local Schur primitive for a
mixed compatible discretization. It eliminates a caller-owned flux block
without changing the flux recovery map:

\[
  Mq+Fu=f,\qquad Gq+Du=g.
\]

The routine returns

\[
  r=M^{-1}f,\qquad X=-M^{-1}F,
  \qquad S=D+GX,\qquad b=g-Gr,
\]

so a condensed solve `S*u=b` recovers the local compatible field as
`q=r+X*u`. FortFEM does not choose a flux space or a global skeleton
numbering: `M`, `F`, `G`, and `D` may come from RT, BDM, Nédélec, mixed IGA,
or another compatible block supplied by the caller.

```fortran
call assemble_compatible_flux_elimination( &
    flux_mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, &
    state_rhs, recovery, recovery_matrix, condensed_matrix, condensed_rhs, &
    status)
call assemble_compatible_flux_elimination_jvp( &
    flux_mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, state_rhs, &
    flux_mass_dot, flux_to_state_dot, state_to_flux_dot, state_matrix_dot, &
    flux_rhs_dot, state_rhs_dot, recovery_dot, recovery_matrix_dot, &
    condensed_matrix_dot, condensed_rhs_dot, status)
call assemble_compatible_flux_elimination_vjp( &
    flux_mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, state_rhs, &
    recovery, recovery_matrix, condensed_matrix, condensed_rhs, recovery_bar, &
    recovery_matrix_bar, condensed_matrix_bar, condensed_rhs_bar, flux_mass_bar, &
    flux_to_state_bar, state_to_flux_bar, state_matrix_bar, flux_rhs_bar, &
    state_rhs_bar, status)
```

The JVP differentiates the implicit mass solves, including changing mass
matrices and both right-hand sides. The VJP uses transpose mass solves and
returns the full real cotangent of every input block. This makes the
operation suitable for fixed-topology optimization and for differentiable
HDG or mixed-IGA assembly while leaving pivot choices and global sparse
ownership to the client.

`test_compatible_flux_elimination` checks the recovery and Schur identities
against direct matrix products, a central-difference JVP, the real VJP dot
identity, and rejection of incompatible output shapes.
