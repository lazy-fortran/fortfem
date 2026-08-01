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
call build_fourier_mode_triad_map( &
    registry, triad_map, missing_count, status)
call build_fourier_mode_padded_registry( &
    registry, padded_registry, status)
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
call apply_fourier_bilinear_product( &
    input_registry, output_registry, coupling, left_values, right_values, &
    product_values, status)
call apply_fourier_bilinear_product_jvp( &
    input_registry, output_registry, coupling, left_values, right_values, &
    coupling_dot, left_values_dot, right_values_dot, product_values_dot, status)
call apply_fourier_bilinear_product_vjp( &
    input_registry, output_registry, coupling, left_values, right_values, &
    product_values_bar, left_values_bar, right_values_bar, coupling_bar, status)
```

`real_packed=.true.` requires every retained `(m,n)` to have its `(-m,-n)`
partner with the same normalization. This is metadata validation only; the
caller still owns the real cosine/sine coefficient packing. Triad closure is
queried explicitly rather than silently padding a truncated mode set.
`build_fourier_mode_triad_map` returns the retained output index for every
ordered input pair and counts omitted sums. A caller can therefore reject an
aliasing-prone set, use `build_fourier_mode_padded_registry` for a deterministic
one-product work set, or deliberately project omitted modes; the neutral
product kernel never makes that policy implicitly. The padded constructor keeps
the original modes first, appends unique ordered pair sums, preserves field
periods, phases, real/conjugate packing, radial powers, and normalizations for
retained modes, and uses `abs(m)` plus unit normalization for newly introduced
modes. It intentionally does not promise closure under repeated products.

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
coefficient symmetry remain caller-owned. `build_fourier_mode_triad_map`
returns the output index for every ordered input pair and an omission count,
so a caller can choose padding, projection, or rejection explicitly. The
one-product padded registry is therefore a work-set constructor, not a hidden
nonlinear truncation rule. Its JVP
differentiates all three factors and its complex VJP uses the same real-part
inner product.

`apply_fourier_bilinear_product` is the de-aliased application form. Its input
arrays use `input_registry`, while the output array uses `output_registry`; the
pair label is looked up only in the output set. Thus a padded work registry can
receive one product without pretending that input modes occupy every padded
slot. Omitted output labels are zero by the explicit projection policy. The
JVP and VJP expose the same distinct-registry map and use the real-part complex
adjoint convention above.

## Independent verification

`test_fourier_mode_registry` checks deep-copy assignment, metadata validation,
conjugate lookup, retained and absent triads, one-product padded closure and
the explicit next-product omission policy, metadata preservation, the
phase/radial analytical formula, central-difference JVP, complex real-part VJP
identity, and malformed real-packed, duplicate, and index inputs.
`test_fourier_vector_product` checks an independent retained-triad map and
contraction,
the full three-factor product-rule JVP and central differences, the complex
real-part VJP identity, and incompatible output-shape rejection.
`test_fourier_bilinear_product` checks the distinct input/output registry
contraction against an independent padded oracle, product-rule JVP and central
differences, the complex real-part VJP identity, and incompatible output-shape
rejection.
