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
call assemble_fourier_vector_product( &
    registry, coupling, left_values, right_values, product_values, status)
call assemble_fourier_vector_product_jvp( &
    registry, coupling, left_values, right_values, coupling_dot, &
    left_values_dot, right_values_dot, product_values_dot, status)
call assemble_fourier_vector_product_vjp( &
    registry, coupling, left_values, right_values, product_values_bar, &
    left_values_bar, right_values_bar, coupling_bar, status)
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

The companion vector product uses a caller-owned coupling tensor and retains
only pairs whose two mode indices sum to another registered mode:

\[
  P_{a,x,k}=\sum_{p+q=k}\sum_{b,c}
      C_{abc}L_{b,x,p}R_{c,x,q}.
\]

It is a pointwise algebraic primitive rather than a model-specific bracket.
Consequently the component counts can represent scalar, vector, tensor,
H(curl), or H(div) coefficients, and a FortSym-generated constitutive block
can supply `C`. Truncated triads are explicit omissions; de-aliasing and
coefficient symmetry remain caller-owned. Its JVP differentiates all three
factors and its complex VJP uses the same real-part inner product.

## Independent verification

`test_fourier_mode_registry` checks deep-copy assignment, metadata validation,
conjugate lookup, retained and absent triads, the phase/radial analytical
formula, central-difference JVP, complex real-part VJP identity, and malformed
real-packed, duplicate, and index inputs.
`test_fourier_vector_product` checks an independent retained-triad contraction,
the full three-factor product-rule JVP and central differences, the complex
real-part VJP identity, and incompatible output-shape rejection.
