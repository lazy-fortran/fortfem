module fortfem_spherical_harmonics
    !! FortFEM adapter for the pinned FortNum spherical-harmonic contract.
    !!
    !! FortNum owns normalization, phase, angular domains, and the pole
    !! convention.  This adapter makes that dependency part of FortFEM's
    !! public numerical API.
    use fortnum_special_spherical, only: &
        spherical_harmonic, spherical_harmonic_theta_derivative, &
        spherical_harmonic_phi_derivative, spherical_harmonic_product_coefficient
    implicit none
    private

    public :: spherical_harmonic
    public :: spherical_harmonic_theta_derivative
    public :: spherical_harmonic_phi_derivative
    public :: spherical_harmonic_product_coefficient

end module fortfem_spherical_harmonics
