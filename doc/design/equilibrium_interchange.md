---
title: Neutral equilibrium interchange contract
---

# Neutral equilibrium interchange contract

`fortfem_equilibrium_interchange` is the boundary between FortFEM's generic
finite-element foundation and an application-owned equilibrium adapter.  It
stores sampled mapped coordinates, physical coordinates, named coefficient and
profile arrays, segmented boundaries, and explicit unit/normalization scales.
The record carries producer and provenance strings, but it does not parse
GEQDSK, COCOS, CHEASE, FreeGS, VMEC, GVEC, DESC, SPEC, GPEC, MARS, or any
other application format.

The schema is deliberately sample based.  A client may populate it from a
license-compatible reader, an analytical manufactured solution, or a common
physical-grid resampler.  The same record can then feed scalar, tensor,
Fourier, IGA, FEM, BEM, or DtN experiments without making those applications
dependencies of FortFEM.

## API

```fortran
type(equilibrium_normalization_t) :: normalization
type(equilibrium_interchange_t) :: data

call initialize_equilibrium_interchange( &
    data, mapped_coordinates, physical_coordinates, coefficient_names, &
    coefficient_values, profile_coordinates, profile_names, profile_values, &
    boundary_offsets, boundary_names, boundary_coordinates, normalization, &
    status)
valid = validate_equilibrium_interchange(data, status)
```

Coordinate arrays have shape `(spatial_dimension, sample_count)` and use the
same active dimension for mapped, physical, and boundary points.  A boundary
is a contiguous point range; `boundary_offsets` has one more entry than
`boundary_names`, starts at one, ends at `size(boundary_coordinates,2)+1`, and
is nondecreasing.  Coefficients are `(coefficient_count, sample_count)` and
profiles are `(profile_count, profile_sample_count)`, with unique nonempty
labels.

Normalization stores unit labels and positive scale factors for length,
magnetic field, pressure, and current.  These are metadata, not hidden unit
conversions.  The application owns the physical interpretation and any
dimensionless normalization convention.

## Differentiation and topology

The record is a fixed-topology interchange object.  Intrinsic assignment makes
a deep copy of all allocatable arrays.  Validation rejects incompatible
shapes, duplicate labels, nonfinite samples, invalid boundary segmentation,
and nonpositive/nonfinite normalization scales before a client assembles a
residual.  Geometry or profile derivatives belong to the external adapter and
can be supplied to FortFEM's residual JVP/VJP contracts separately.  Changing
the number or ordering of samples is a topology/schema event, not a silently
differentiated perturbation.

The manufactured test uses toroidal coordinates

\[
 x=(R+\rho\cos\theta)\cos\zeta,\quad
 y=(R+\rho\cos\theta)\sin\zeta,\quad
 z=\rho\sin\theta,
\]

and independently checks coordinate, coefficient, profile, boundary,
normalization, copy, and rejection behavior.  This keeps the contract
license-safe while providing the common-grid oracle needed by future CHEASE
and FreeGS adapters.
