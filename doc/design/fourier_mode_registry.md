---
title: Fourier mode registry
---

# Fourier mode registry

`fourier_mode_registry_t` is the fixed-topology metadata contract shared by
Fourier-FEM, multipatch IGA, and caller-defined toroidal operators. It keeps
mode numbering and normalization out of individual element kernels.

For a retained pair `(m,n)`, the registry uses the documented phase
convention

\[
  \chi_{mn}(r,\theta,\phi)=a_{mn}r^{p_{mn}}
    \exp\!\left(i\left[m(\theta+\theta_0)+nN_{\rm fp}
    (\phi+\phi_0)\right]\right),
\]

where `N_fp` is the positive number of field periods, `a_mn` is the caller's
positive normalization, and `p_mn` is the caller-selected nonnegative radial
regularity power. The radial factor is evaluated with its analytic derivative
at the axis (`p=0` and `p=1` have their regular one-sided values).

## API

```fortran
call initialize_fourier_mode_registry( &
    registry, poloidal_modes, toroidal_modes, field_periods, &
    poloidal_phase, toroidal_phase, real_packed, radial_powers, &
    normalization, status)
valid = validate_fourier_mode_registry(registry, status)
index = find_fourier_mode(registry, poloidal_mode, toroidal_mode)
conjugate = fourier_mode_conjugate_index(registry, index, status)
closed = fourier_mode_triad_closed(registry, first, second, output, status)
call evaluate_fourier_mode( &
    registry, index, radius, theta, phi, value, radial_derivative, &
    theta_derivative, phi_derivative, status)
call evaluate_fourier_mode_jvp( &
    registry, index, radius, theta, phi, radius_dot, theta_dot, phi_dot, &
    value_dot, status)
call evaluate_fourier_mode_vjp( &
    registry, index, radius, theta, phi, value_bar, radius_bar, theta_bar, &
    phi_bar, status)
```

`real_packed=.true.` requires every retained `(m,n)` to have its `(-m,-n)`
partner with the same normalization. This is metadata validation only; the
caller still owns the real cosine/sine coefficient packing. Triad closure is
queried explicitly rather than silently padding a truncated mode set.

The value, fixed-topology JVP, and real-coordinate VJP use the complex
real-part convention

\[
  \langle z_{\rm bar},\dot z\rangle=\operatorname{Re}
  (\overline{z_{\rm bar}}\,\dot z).
\]

No equilibrium profiles, COCOS/GEQDSK conventions, or plasma closure are
encoded here. A future torus-harmonic or FortNum special-function adapter can
use the registry's indices and replace only the radial factor while retaining
its phase, packing, and derivative contracts.

## Independent verification

`test_fourier_mode_registry` checks deep-copy assignment, metadata validation,
conjugate lookup, retained and absent triads, the phase/radial analytical
formula, central-difference JVP, complex real-part VJP identity, and malformed
real-packed, duplicate, and index inputs.
