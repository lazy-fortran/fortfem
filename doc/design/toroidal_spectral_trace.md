---
title: Toroidal spectral trace contract
---

# Toroidal spectral trace contract

`evaluate_toroidal_spectral_trace` is a neutral modal boundary operator for a
constant-η toroidal surface.  It is a building block for a NESTOR-like
periodic Green operator; it does not select an equilibrium, coil model, field
period, or application normalization.

For a mode `(degree_index, order)`, FortNum supplies the Hobson-normalized
half-integer toroidal function (R=P_{n-1/2}^m) or (R=Q_{n-1/2}^m).  The
complex basis used by the adapter is

\[
  \psi_{nm}(\eta,\theta,\phi)
  =\sqrt{\cosh\eta-\cos\theta}\,R(\cosh\eta)
    \exp\{i(n\theta+m\phi)\}.
\]

The scale is the focal length (a>0).  The returned normal derivative uses
the outward normal of the region η>η₀, namely the direction of decreasing
η:

\[
  \partial_n\psi
  =-\frac{\cosh\eta-\cos\theta}{a}\,\partial_\eta\psi.
\]

The modal coefficient, scale, and angular-coordinate JVP/VJP actions are
fixed-topology analytical products.  The test compares the P branch against
the existing single toroidal harmonic and DtN oracle, and checks central
reassembly and the real-part complex adjoint identity.  Q modes require the
valid FortNum degree/order region; zero-mode handling, periodic sampling,
convolution, and truncation estimates belong to the higher-level spectral
operator and remain caller policy.
