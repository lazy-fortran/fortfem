---
title: Axis-regular Fourier mode metadata
---

# Axis-regular Fourier mode metadata

`fortfem_axis_regular_fourier_modes` is a neutral metadata contract for a
scalar Fourier coefficient on a polar or toroidal magnetic axis.  For a
poloidal label (m), a smooth scalar coefficient has radial powers

\[
  p = |m|, |m|+2, |m|+4, \ldots .
\]

The minimum power is therefore `abs(m)` and the required parity is
`modulo(abs(m),2)`.  `axis_regular_mode_requirements` reports those two values
and returns an invalid status when the supplied power is negative, below the
minimum, or has the wrong parity.  This is a scalar basis rule; vector and
tensor components may use shifted effective labels and must perform that
choice in their own space metadata.

## Registry table

```fortran
call build_axis_regular_mode_table(registry, table, status)
valid = validate_axis_regular_mode_table(table, status)
```

The builder consumes the existing `fourier_mode_registry_t`, retains field
period and real-packing metadata, and returns one record per retained mode.
Each record reports the original registry index, `(m,n)`, supplied radial
power, minimum power, parity, axis-regular flag, and both registry and ordered
conjugate indices.  A real-packed registry must contain every conjugate (as
required by the Fourier registry itself) and conjugate radial powers must
match.

The output order is independent of input array order.  A canonical
representative of each conjugate pair is placed first, followed immediately by
its conjugate; the key is `(m,n)` for `m>0` (or `m=0,n>=0`), otherwise the
negated pair, followed by orientation.  This keeps complex packing and
field-period mode lookup deterministic without selecting a plasma model.

The table only validates fixed-topology metadata.  It does not choose profiles,
enforce force balance, or import DESC, VMEC, or any plasma-specific format.
`evaluate_axis_regular_radial_basis` applies the same rule to a caller-selected
finite polynomial and supplies coefficient/radial-coordinate derivatives;
`nested_surface_geometry` can consume the common mode metadata independently.

## Independent verification

`test_axis_regular_fourier_modes` checks the scalar (m=0) and (m>0)
contracts, rejection of parity and minimum-power violations, deterministic
ordering under a shuffled registry, conjugate adjacency for real-packed
modes, and the intentional missing-conjugate policy for a complex registry.
