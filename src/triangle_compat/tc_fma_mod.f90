module tc_fma_mod
    !! Portable correctly rounded fused multiply-add.  Fortran 2023
    !! ieee_fma is missing from gfortran < 13, so bind to the C99 libm fma,
    !! which is exact (a single rounding) on every conforming libm.
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: tc_fma

    interface
        pure function c_fma(x, y, z) bind(C, name="fma") result(r)
            import :: c_double
            real(c_double), value, intent(in) :: x, y, z
            real(c_double) :: r
        end function c_fma
    end interface

contains

    pure elemental function tc_fma(a, b, c) result(r)
        real(dp), intent(in) :: a, b, c
        real(dp) :: r

        r = real(c_fma(real(a, c_double), real(b, c_double), real(c, c_double)), dp)
    end function tc_fma

end module tc_fma_mod
