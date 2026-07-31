module fortfem_toroidal_harmonics
    !! FortFEM adapter for the pinned FortNum toroidal P/Q contract.
    !!
    !! FortNum owns the Hobson normalization and half-integer recurrence.
    !! This module keeps those analytical torus-harmonic functions available
    !! through the FortFEM public API without selecting an application law.
    use fortnum_special_toroidal, only: &
        toroidal_p, toroidal_q, toroidal_p_derivative, toroidal_q_derivative
    implicit none
    private

    public :: toroidal_p
    public :: toroidal_q
    public :: toroidal_p_derivative
    public :: toroidal_q_derivative

end module fortfem_toroidal_harmonics
