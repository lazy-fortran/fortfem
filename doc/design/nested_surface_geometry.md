% Nested-surface geometry contract
% FortFEM

# Scope

`fortfem_nested_surface_geometry` is a small geometry-only building block for
external nested-surface equilibrium and response clients.  It does not read
EQDSK/VMEC/DESC files, choose a pressure or current model, enforce a flux
normalization, or implement a plasma equation.

The caller supplies a `fourier_mode_registry_t` and a complex coefficient array
with shape `(3, nmode)`.  The three rows are the real fields

\[
  (R,Z,\lambda)_c(\rho,\theta,\zeta)
  = \operatorname{Re}\sum_k c_{c k}
      b_k(\rho,\theta,\zeta),
\]

where `b_k` is the existing FortFEM Fourier/radial registry basis, including
its field-period, phase, normalization, and radial-power metadata.  This
keeps axis regularity and mode ordering in one shared registry instead of
duplicating a Fourier implementation.

The inverse-coordinate map uses the neutral convention

\[
  \phi=\zeta+\lambda,
  \qquad
  (x,y,z)=(R\cos\phi,R\sin\phi,Z).
\]

The convention is intentionally explicit: an adapter using a different
physical toroidal-angle normalization must transform its samples or provide a
different geometry layer.

# Returned diagnostics

For every fixed sample `(rho,theta,zeta)`, the evaluator returns:

* inverse coordinates `(R,Z,lambda)`;
* Cartesian coordinates;
* the inverse-coordinate Jacobian
  `d(R,Z,lambda)/d(rho,theta,zeta)`;
* the Cartesian Jacobian;
* the covariant metric `J^T J`;
* the signed Cartesian volume Jacobian `det(J)`.

Zero volume at a coordinate singularity is reported as a diagnostic; it is
not silently replaced by a positive value.  The routine only rejects
non-finite data and negative radial coordinates.  A client remains responsible
for deciding whether its selected axis basis makes a nonsingular physical
surface.

# Differentiation

The JVP and VJP differentiate coefficients at fixed sample coordinates and
fixed mode topology.  The VJP includes all returned quantities (coordinates,
both Jacobians, metric, and volume determinant).  For a complex coefficient
increment `dc`, the real pairing is

\[
  \langle \bar c,dc\rangle=\operatorname{Re}\sum_k
  \overline{\bar c_k}\,dc_k.
\]

Mode activation, radial powers, field-period metadata, and phase conventions
are discrete topology/metadata choices.  Changing them requires rebuilding
the registry and is outside the fixed-topology derivative contract.

# Provenance and intended composition

The radial Fourier basis and phase convention are delegated to
`fortfem_fourier_mode_registry`; its analytical radial and angular derivatives
are therefore the single source of truth.  The Cartesian toroidal embedding is
the same `phi=zeta+lambda` convention exposed by the neutral geometry layer,
not a claim about any one equilibrium code.  External VMEC/GVEC/DESC-like
optimizers can use the returned geometry and derivatives to build their own
force residuals, objectives, flux constraints, Boozer diagnostics, or
free-boundary coupling.  The implementation is also usable by fitted IGA,
Fourier-FEM, BEM, DtN, and virtual-casing clients that need physical sample
coordinates without importing plasma-specific formats.

The independent test covers a manufactured modal oracle, field-period
periodicity, the radial-power axis contract, central-difference JVP checks, and
the coefficient VJP dot-product identity.
