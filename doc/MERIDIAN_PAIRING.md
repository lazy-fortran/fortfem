# Meridian 1-form/2-form pairing

`meridian_mass_weighted` accepts optional `trial_space`. Its test space remains
`space`, and the resulting matrix has shape `(sp(space)%ndof,
sp(trial_space)%ndof)`. The default is the original square mass. Both physical
basis vectors, their orientation signs, the cylindrical measure and an optional
scalar weight are included. The contraction and its JVP/VJP are retained in
`tools/codegen/app/gen_triangle_piola_products.f90` and checked-in output.
Element triplets are assembled with OpenMP and no shared scatter.

`meridian_fourier_pairing(fs,modes,m,status)` returns the rectangular real
Fourier V1-V2 pairing. The mean block has weight one and cosine/sine blocks
have weight one half. It is metric invariant under paired covariant 1-form and
contravariant 2-form Piola transformations. For a discrete curl C, `P C` is the
helicity form, symmetric on variations with zero tangential trace. No global
symmetry is asserted for unrestricted boundary traces.

`meridian_fourier_diagonal_mass(fs,space,modes,m,status)` lifts an existing
meridional mass using the same exact Fourier weights. It avoids storing
roundoff cross-mode couplings from trigonometric quadrature in an identity
metric. General phi-dependent tensors retain `meridian_fourier_mass`.

The focused `test_meridian_fourier_curl` compares the mixed pairing on smooth
polynomial-trigonometric fields with independent tensor-product quadrature of
analytic A and curl A. It also compares the diagonal mass with independent
uniform-phi integration and retains the original curl/gradient checks.
