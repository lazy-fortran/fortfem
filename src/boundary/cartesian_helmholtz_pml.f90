module fortfem_cartesian_helmholtz_pml
    !! Complex-coordinate PML tensors for scalar Helmholtz and Maxwell.
    !!
    !! With diagonal stretch S=diag(sx,sy,sz), the H1 weak-form tensor is
    !! det(S) S^-1 S^-T.  The H(curl) curl tensor is its inverse, while the
    !! Maxwell mass tensor is the H1 tensor.  This is the exact-sequence/Piola
    !! construction of Matuszyk & Demkowicz, Computational Mechanics 51 (2013),
    !! DOI:10.1007/s00466-012-0702-1.
    use fortfem_kinds, only: dp
    use fortfem_generated_cartesian_curl_curl_pml_jvp, only: &
        generated_cartesian_curl_curl_pml_jvp
    use fortfem_generated_cartesian_curl_curl_pml_vjp, only: &
        generated_cartesian_curl_curl_pml_vjp
    use fortfem_generated_cartesian_scalar_pml_jvp, only: &
        generated_cartesian_scalar_pml_jvp
    use fortfem_generated_cartesian_scalar_pml_vjp, only: &
        generated_cartesian_scalar_pml_vjp
    implicit none
    private

    public :: cartesian_curl_curl_pml_coefficients
    public :: cartesian_curl_curl_pml_coefficients_jvp
    public :: cartesian_curl_curl_pml_coefficients_vjp
    public :: cartesian_scalar_helmholtz_pml_coefficients
    public :: cartesian_scalar_helmholtz_pml_coefficients_jvp
    public :: cartesian_scalar_helmholtz_pml_coefficients_vjp

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

    pure subroutine cartesian_scalar_helmholtz_pml_coefficients_jvp( &
            stretch, stretch_dot, gradient_coefficient_dot, &
            mass_coefficient_dot, status)
        complex(dp), intent(in) :: stretch(3), stretch_dot(3)
        complex(dp), intent(out) :: gradient_coefficient_dot(3)
        complex(dp), intent(out) :: mass_coefficient_dot
        integer, intent(out) :: status

        gradient_coefficient_dot = cmplx(0.0_dp, 0.0_dp, dp)
        mass_coefficient_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (any(abs(stretch) <= tiny(1.0_dp))) return
        call generated_cartesian_scalar_pml_jvp( &
            stretch(1), stretch(2), stretch(3), stretch_dot(1), &
            stretch_dot(2), stretch_dot(3), gradient_coefficient_dot(1), &
            gradient_coefficient_dot(2), gradient_coefficient_dot(3), &
            mass_coefficient_dot)
        status = 0
    end subroutine cartesian_scalar_helmholtz_pml_coefficients_jvp

    pure subroutine cartesian_scalar_helmholtz_pml_coefficients_vjp( &
            stretch, gradient_coefficient_bar, mass_coefficient_bar, &
            stretch_bar, status)
        complex(dp), intent(in) :: stretch(3), gradient_coefficient_bar(3)
        complex(dp), intent(in) :: mass_coefficient_bar
        complex(dp), intent(out) :: stretch_bar(3)
        integer, intent(out) :: status

        complex(dp) :: holomorphic_bar(3)

        stretch_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (any(abs(stretch) <= tiny(1.0_dp))) return
        call generated_cartesian_scalar_pml_vjp( &
            stretch(1), stretch(2), stretch(3), &
            conjg(gradient_coefficient_bar(1)), &
            conjg(gradient_coefficient_bar(2)), &
            conjg(gradient_coefficient_bar(3)), conjg(mass_coefficient_bar), &
            holomorphic_bar(1), holomorphic_bar(2), holomorphic_bar(3))
        stretch_bar = conjg(holomorphic_bar)
        status = 0
    end subroutine cartesian_scalar_helmholtz_pml_coefficients_vjp

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

    pure subroutine cartesian_curl_curl_pml_coefficients_jvp( &
            stretch, stretch_dot, curl_coefficient_dot, &
            mass_coefficient_dot, status)
        complex(dp), intent(in) :: stretch(3), stretch_dot(3)
        complex(dp), intent(out) :: curl_coefficient_dot(3)
        complex(dp), intent(out) :: mass_coefficient_dot(3)
        integer, intent(out) :: status

        curl_coefficient_dot = cmplx(0.0_dp, 0.0_dp, dp)
        mass_coefficient_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (any(abs(stretch) <= tiny(1.0_dp))) return
        call generated_cartesian_curl_curl_pml_jvp( &
            stretch(1), stretch(2), stretch(3), stretch_dot(1), &
            stretch_dot(2), stretch_dot(3), curl_coefficient_dot(1), &
            curl_coefficient_dot(2), curl_coefficient_dot(3), &
            mass_coefficient_dot(1), mass_coefficient_dot(2), &
            mass_coefficient_dot(3))
        status = 0
    end subroutine cartesian_curl_curl_pml_coefficients_jvp

    pure subroutine cartesian_curl_curl_pml_coefficients_vjp( &
            stretch, curl_coefficient_bar, mass_coefficient_bar, &
            stretch_bar, status)
        complex(dp), intent(in) :: stretch(3), curl_coefficient_bar(3)
        complex(dp), intent(in) :: mass_coefficient_bar(3)
        complex(dp), intent(out) :: stretch_bar(3)
        integer, intent(out) :: status

        complex(dp) :: holomorphic_bar(3)

        stretch_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (any(abs(stretch) <= tiny(1.0_dp))) return
        call generated_cartesian_curl_curl_pml_vjp( &
            stretch(1), stretch(2), stretch(3), &
            conjg(curl_coefficient_bar(1)), conjg(curl_coefficient_bar(2)), &
            conjg(curl_coefficient_bar(3)), conjg(mass_coefficient_bar(1)), &
            conjg(mass_coefficient_bar(2)), conjg(mass_coefficient_bar(3)), &
            holomorphic_bar(1), holomorphic_bar(2), holomorphic_bar(3))
        stretch_bar = conjg(holomorphic_bar)
        status = 0
    end subroutine cartesian_curl_curl_pml_coefficients_vjp

end module fortfem_cartesian_helmholtz_pml
