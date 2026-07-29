module fortfem_cartesian_helmholtz_pml
    !! Complex-coordinate PML tensors for scalar Helmholtz and Maxwell.
    !!
    !! With diagonal stretch S=diag(sx,sy,sz), the H1 weak-form tensor is
    !! det(S) S^-1 S^-T.  The H(curl) curl tensor is its inverse, while the
    !! Maxwell mass tensor is the H1 tensor.  This is the exact-sequence/Piola
    !! construction of Matuszyk & Demkowicz, Computational Mechanics 51 (2013),
    !! DOI:10.1007/s00466-012-0702-1.
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: cartesian_curl_curl_pml_coefficients
    public :: cartesian_scalar_helmholtz_pml_coefficients

contains

    pure subroutine cartesian_scalar_helmholtz_pml_coefficients( &
            stretch, gradient_coefficient, mass_coefficient, status)
        complex(dp), intent(in) :: stretch(3)
        complex(dp), intent(out) :: gradient_coefficient(3)
        complex(dp), intent(out) :: mass_coefficient
        integer, intent(out) :: status

        gradient_coefficient = cmplx(0.0_dp, 0.0_dp, dp)
        mass_coefficient = cmplx(0.0_dp, 0.0_dp, dp)
        status = 0
        if (any(abs(stretch) <= tiny(1.0_dp))) then
            status = 1
            return
        end if

        gradient_coefficient(1) = stretch(2)*stretch(3)/stretch(1)
        gradient_coefficient(2) = stretch(1)*stretch(3)/stretch(2)
        gradient_coefficient(3) = stretch(1)*stretch(2)/stretch(3)
        mass_coefficient = product(stretch)
    end subroutine cartesian_scalar_helmholtz_pml_coefficients

    pure subroutine cartesian_curl_curl_pml_coefficients( &
            stretch, curl_coefficient, mass_coefficient, status)
        complex(dp), intent(in) :: stretch(3)
        complex(dp), intent(out) :: curl_coefficient(3)
        complex(dp), intent(out) :: mass_coefficient(3)
        integer, intent(out) :: status

        complex(dp) :: scalar_mass

        call cartesian_scalar_helmholtz_pml_coefficients( &
            stretch, mass_coefficient, scalar_mass, status)
        curl_coefficient = cmplx(0.0_dp, 0.0_dp, dp)
        if (status /= 0) return
        curl_coefficient = 1.0_dp/mass_coefficient
    end subroutine cartesian_curl_curl_pml_coefficients

end module fortfem_cartesian_helmholtz_pml
