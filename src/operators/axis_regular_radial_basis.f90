module fortfem_axis_regular_radial_basis
    !! Differentiable finite radial polynomials for one scalar Fourier mode.
    !!
    !! The caller selects distinct, increasing powers and complex coefficients.
    !! Every power must satisfy the same scalar axis contract used by
    !! `build_axis_regular_mode_table`: p >= abs(m), with matching parity.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_axis_regular_fourier_modes, only: &
        axis_regular_mode_requirements
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_axis_regular_radial_basis
    public :: evaluate_axis_regular_radial_basis_jvp
    public :: evaluate_axis_regular_radial_basis_vjp

contains

    subroutine evaluate_axis_regular_radial_basis( &
            poloidal_mode, radial_powers, coefficients, rho, values, status)
        integer, intent(in) :: poloidal_mode
        integer, intent(in) :: radial_powers(:)
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: rho(:)
        complex(dp), intent(out) :: values(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: term

        values = cmplx(0.0_dp, 0.0_dp, dp)
        call validate_radial_basis_inputs( &
            poloidal_mode, radial_powers, coefficients, rho, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(values) /= size(rho)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "axis-regular radial basis has an incompatible output shape")
            return
        end if

        do term = 1, size(radial_powers)
            values = values + coefficients(term)*rho**radial_powers(term)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_axis_regular_radial_basis

    subroutine evaluate_axis_regular_radial_basis_jvp( &
            poloidal_mode, radial_powers, coefficients, rho, coefficients_dot, &
            rho_dot, values_dot, status)
        !! Fixed-selector tangent with respect to coefficients and samples.
        integer, intent(in) :: poloidal_mode
        integer, intent(in) :: radial_powers(:)
        complex(dp), intent(in) :: coefficients(:), coefficients_dot(:)
        real(dp), intent(in) :: rho(:), rho_dot(:)
        complex(dp), intent(out) :: values_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: term
        real(dp) :: radial_derivative(size(rho))

        values_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call validate_radial_basis_inputs( &
            poloidal_mode, radial_powers, coefficients, rho, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(coefficients_dot) /= size(coefficients) .or. &
                size(rho_dot) /= size(rho) .or. &
                size(values_dot) /= size(rho)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "axis-regular radial-basis JVP has incompatible shapes")
            return
        end if
        if (.not. finite_complex(coefficients_dot) .or. &
                .not. finite_real(rho_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "axis-regular radial-basis JVP has non-finite tangents")
            return
        end if

        do term = 1, size(radial_powers)
            call evaluate_power_derivative( &
                radial_powers(term), rho, radial_derivative)
            values_dot = values_dot + &
                coefficients_dot(term)*rho**radial_powers(term) + &
                coefficients(term)*radial_derivative*rho_dot
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_axis_regular_radial_basis_jvp

    subroutine evaluate_axis_regular_radial_basis_vjp( &
            poloidal_mode, radial_powers, coefficients, rho, values_bar, &
            coefficients_bar, rho_bar, status)
        !! Reverse action under complex real-part and real Euclidean pairings.
        integer, intent(in) :: poloidal_mode
        integer, intent(in) :: radial_powers(:)
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: rho(:)
        complex(dp), intent(in) :: values_bar(:)
        complex(dp), intent(out) :: coefficients_bar(:)
        real(dp), intent(out) :: rho_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: term
        real(dp) :: radial_derivative(size(rho))

        coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        rho_bar = 0.0_dp
        call validate_radial_basis_inputs( &
            poloidal_mode, radial_powers, coefficients, rho, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(values_bar) /= size(rho) .or. &
                size(coefficients_bar) /= size(coefficients) .or. &
                size(rho_bar) /= size(rho)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "axis-regular radial-basis VJP has incompatible shapes")
            return
        end if
        if (.not. finite_complex(values_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "axis-regular radial-basis VJP has non-finite cotangents")
            return
        end if

        do term = 1, size(radial_powers)
            coefficients_bar(term) = &
                sum(values_bar*rho**radial_powers(term))
            call evaluate_power_derivative( &
                radial_powers(term), rho, radial_derivative)
            rho_bar = rho_bar + real( &
                conjg(values_bar)*coefficients(term), dp)*radial_derivative
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_axis_regular_radial_basis_vjp

    subroutine validate_radial_basis_inputs( &
            poloidal_mode, radial_powers, coefficients, rho, status)
        integer, intent(in) :: poloidal_mode
        integer, intent(in) :: radial_powers(:)
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: rho(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: term, minimum_power, parity
        logical :: regular

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "axis-regular radial basis has invalid inputs")
        if (size(radial_powers) < 1 .or. size(rho) < 1) return
        if (size(coefficients) /= size(radial_powers)) return
        if (.not. finite_complex(coefficients)) return
        if (.not. finite_real(rho)) return
        if (any(rho < 0.0_dp)) return
        do term = 2, size(radial_powers)
            if (radial_powers(term) <= radial_powers(term - 1)) return
        end do
        do term = 1, size(radial_powers)
            call axis_regular_mode_requirements( &
                poloidal_mode, radial_powers(term), minimum_power, parity, &
                regular, status)
            if (.not. regular) return
            if (status%code /= FORTSPARSE_OK) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_radial_basis_inputs

    subroutine evaluate_power_derivative(power, rho, derivative)
        integer, intent(in) :: power
        real(dp), intent(in) :: rho(:)
        real(dp), intent(out) :: derivative(:)

        if (power == 0) then
            derivative = 0.0_dp
        else
            derivative = real(power, dp)*rho**(power - 1)
        end if
    end subroutine evaluate_power_derivative

    logical function finite_real(values) result(finite)
        real(dp), intent(in) :: values(:)

        finite = all(ieee_is_finite(values))
    end function finite_real

    logical function finite_complex(values) result(finite)
        complex(dp), intent(in) :: values(:)

        finite = all(ieee_is_finite(real(values, dp)))
        if (.not. finite) return
        finite = all(ieee_is_finite(aimag(values)))
    end function finite_complex

end module fortfem_axis_regular_radial_basis
