# Physical harmonic field sectors

[V] `fortfem_harmonic_field_energy` couples caller-supplied physical harmonic
representatives to a scalar energy form. For a genus-g handlebody there are g
volume sectors on each side of the boundary; the 2g surface cohomology basis
alone does not construct these fields. The caller must supply curl-free,
divergence-free representatives on the actual domain, with their physical
period or cut-flux functionals and energy integrals. Exterior representatives
must also have finite energy at infinity.

## Assembly and prescribed coordinates

[V] `assemble_field_gram(weights,fields,gram,status)` integrates the Euclidean
field inner products with nonnegative volume quadrature weights. Storage is
`fields(3,npoints,nfields)`, with no allocations in the quadrature loop.
The pointwise product and its JVP/VJP are emitted by FortSym. The quadrature
reduction uses OpenMP even when there is only one harmonic field.

[V] Partition the resulting Gram matrix into scalar, cross and harmonic
blocks S, C and G. Alternatively obtain S from a scalar DtN operator, C from
the corresponding boundary integral and G from an independent harmonic-field
construction. All three blocks must describe the same domain and field
representatives. For scalar field `grad(u)`, integration by parts gives
`C_ij = integral_boundary u_i h_j.n_domain`. On an exterior domain the normal
at the inner boundary is opposite to FortFEM's outward-from-interior normal,
so the exterior cross block uses a minus sign. Quadrature and field
construction errors remain the caller's responsibility.

```fortran
call couple_harmonic_sector(S, C, G, coordinates, prescribed, u, &
    coefficients, load, energy, status)
```

[V] The square g-by-g coordinate map contains the chosen independent physical
functionals applied to each representative. Rows may prescribe periods or
cut fluxes, with their original units. The routine solves
`coordinates*coefficients = prescribed` and returns
`energy = u^T S u + 2 u^T C coefficients + coefficients^T G coefficients`.
The scalar load is `S*u + C*coefficients`, half the energy derivative with
respect to u at fixed prescribed coordinates. There is no factor of one half
in the returned energy; permeability factors belong to the caller.

[V] Prescribing both a period and its dependent cut flux as independent data
overdetermines a sector. Choose g independent coordinates. Flux from the
scalar field must be subtracted from the requested total flux before calling
this routine. Coordinate functionals are held fixed during evaluation; this
API does not differentiate a moving cycle or cut.

[V] Status 1 rejects invalid shapes, nonfinite inputs, negative quadrature
weights, nonsymmetric energy blocks or nonpositive harmonic diagonal entries.
Status 2 rejects a singular coordinate solve, nonfinite output or negative
evaluated total energy beyond roundoff. Failed coupling outputs are zero.
These guards do not certify positive semidefiniteness of arbitrary supplied
blocks. For quadrature-assembled blocks it follows from the nonnegative
weights; externally assembled blocks need independent consistency checks.

## Exact torus oracle and limits

[V] On the solid torus `(R-R0)^2+z^2<a^2`, with `R0>a`, the field
`h=grad(phi)=(-y,x,0)/(x*x+y*y)` is tangent, curl-free and divergence-free.
Its toroidal period is `2*pi`; its meridional cut flux is
`F=2*pi*(R0-sqrt(R0^2-a^2))`. Integrating `|h|^2=1/R^2` over the torus gives
`E_h=2*pi*F`. One derivation integrates phi first, then uses
`integral_0^(2*pi) dtheta/(R0+r*cos(theta)) = 2*pi/sqrt(R0^2-r^2)`
and integrates `r dr` from zero to a. FortSym checks the Cartesian divergence,
curl and energy-density identities symbolically, away from the z axis.

[V] The independent native test integrates the volume Gram matrix for
`grad(x)` and `H=h+grad(x)`. With torus volume V, the exact blocks are
`S=V`, `C=V`, `G=E_h+V`. This deliberately non-tangent representative makes
omitting or changing the sign of the cross term observable. Prescribing
amplitude c and scalar coefficient `u=-c` gives `c^2*E_h` and zero scalar
load. Separate tests prescribe circulation and meridional flux, and exercise
a coupled two-sector coordinate map and generated derivative leaves.

[V] This is a physical energy/coordinate coupling primitive with a bounded
solid-torus oracle. It does not construct exterior harmonic fields for an
arbitrary surface, establish their decay, or supply a general poloidal
exterior oracle. The nondecaying axial line-current field used on the bounded
torus must not be used as an exterior finite-energy representative.
