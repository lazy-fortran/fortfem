module fortfem_fourier_mode_expansion
    !! Neutral coefficient-to-point Fourier/radial synthesis.
    !!
    !! A registry supplies the retained labels, radial powers, phases, and
    !! normalizations.  This layer contracts caller-owned complex mode
    !! coefficients with all registry basis functions at one coordinate:
    !!
    !!   u(rho,theta,phi) = sum_m c_m b_m(rho,theta,phi).
    !!
    !! It is deliberately independent of a geometry, field interpretation, or
    !! equilibrium model.  The same contraction is useful for scalar, vector,
    !! and tensor blocks when called component by component.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_fourier_mode_registry, only: &
        evaluate_fourier_mode, fourier_mode_registry_t, &
        validate_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_fourier_mode_expansion
    public :: evaluate_fourier_mode_expansion_jvp
    public :: evaluate_fourier_mode_expansion_vjp
    public :: evaluate_fourier_mode_expansion_hessian
    public :: evaluate_fourier_mode_expansion_hvp

contains

    subroutine evaluate_fourier_mode_expansion( &
            registry, coefficients, radius, theta, phi, value, &
            radial_derivative, theta_derivative, phi_derivative, status)
        !! Evaluate all retained registry modes at one cylindrical point.
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: radius, theta, phi
        complex(dp), intent(out) :: value, radial_derivative
        complex(dp), intent(out) :: theta_derivative, phi_derivative
        type(fortsparse_status_t), intent(out) :: status

        integer :: mode
        complex(dp) :: mode_value, mode_radial, mode_theta, mode_phi

        value = cmplx(0.0_dp, 0.0_dp, dp)
        radial_derivative = value
        theta_derivative = value
        phi_derivative = value
        if (.not. validate_expansion_inputs( &
                registry, coefficients, radius, theta, phi, status)) return
        do mode = 1, size(coefficients)
            call evaluate_fourier_mode( &
                registry, mode, radius, theta, phi, mode_value, mode_radial, &
                mode_theta, mode_phi, status)
            if (status%code /= FORTSPARSE_OK) return
            value = value + coefficients(mode)*mode_value
            radial_derivative = radial_derivative + coefficients(mode)*mode_radial
            theta_derivative = theta_derivative + coefficients(mode)*mode_theta
            phi_derivative = phi_derivative + coefficients(mode)*mode_phi
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fourier_mode_expansion

    subroutine evaluate_fourier_mode_expansion_jvp( &
            registry, coefficients, radius, theta, phi, coefficients_dot, &
            radius_dot, theta_dot, phi_dot, value_dot, status)
        !! Apply the exact coefficient/coordinate tangent of the expansion.
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: radius, theta, phi
        complex(dp), intent(in) :: coefficients_dot(:)
        real(dp), intent(in) :: radius_dot, theta_dot, phi_dot
        complex(dp), intent(out) :: value_dot
        type(fortsparse_status_t), intent(out) :: status

        integer :: mode
        complex(dp) :: mode_value, mode_radial, mode_theta, mode_phi

        value_dot = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. validate_expansion_inputs( &
                registry, coefficients, radius, theta, phi, status)) return
        if (size(coefficients_dot) /= size(coefficients) .or. &
                .not. finite_complex(coefficients_dot) .or. &
                .not. ieee_is_finite(radius_dot) .or. &
                .not. ieee_is_finite(theta_dot) .or. &
                .not. ieee_is_finite(phi_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier mode expansion JVP has invalid increments")
            return
        end if
        do mode = 1, size(coefficients)
            call evaluate_fourier_mode( &
                registry, mode, radius, theta, phi, mode_value, mode_radial, &
                mode_theta, mode_phi, status)
            if (status%code /= FORTSPARSE_OK) return
            value_dot = value_dot + coefficients_dot(mode)*mode_value + &
                coefficients(mode)*(mode_radial*radius_dot + &
                mode_theta*theta_dot + mode_phi*phi_dot)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fourier_mode_expansion_jvp

    subroutine evaluate_fourier_mode_expansion_vjp( &
            registry, coefficients, radius, theta, phi, value_bar, &
            coefficients_bar, radius_bar, theta_bar, phi_bar, status)
        !! Apply the real-part complex adjoint of the expansion.
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: radius, theta, phi
        complex(dp), intent(in) :: value_bar
        complex(dp), intent(out) :: coefficients_bar(:)
        real(dp), intent(out) :: radius_bar, theta_bar, phi_bar
        type(fortsparse_status_t), intent(out) :: status

        integer :: mode
        complex(dp) :: mode_value, mode_radial, mode_theta, mode_phi

        coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        radius_bar = 0.0_dp
        theta_bar = 0.0_dp
        phi_bar = 0.0_dp
        if (.not. validate_expansion_inputs( &
                registry, coefficients, radius, theta, phi, status)) return
        if (size(coefficients_bar) /= size(coefficients) .or. &
                .not. finite_complex([value_bar])) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier mode expansion VJP has invalid cotangent shape/value")
            return
        end if
        do mode = 1, size(coefficients)
            call evaluate_fourier_mode( &
                registry, mode, radius, theta, phi, mode_value, mode_radial, &
                mode_theta, mode_phi, status)
            if (status%code /= FORTSPARSE_OK) return
            coefficients_bar(mode) = value_bar*conjg(mode_value)
            radius_bar = radius_bar + real( &
                conjg(value_bar)*coefficients(mode)*mode_radial, dp)
            theta_bar = theta_bar + real( &
                conjg(value_bar)*coefficients(mode)*mode_theta, dp)
            phi_bar = phi_bar + real( &
                conjg(value_bar)*coefficients(mode)*mode_phi, dp)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fourier_mode_expansion_vjp

    subroutine evaluate_fourier_mode_expansion_hessian( &
            registry, coefficients, radius, theta, phi, radial_second, &
            radial_theta, radial_phi, theta_second, theta_phi, phi_second, &
            status)
        !! Evaluate the six symmetric coordinate Hessian components.
        !!
        !! The fixed-topology coordinates are ordered as (rho, theta, phi).
        !! No coefficient or geometry convention beyond the registry is used.
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: radius, theta, phi
        complex(dp), intent(out) :: radial_second, radial_theta, radial_phi
        complex(dp), intent(out) :: theta_second, theta_phi, phi_second
        type(fortsparse_status_t), intent(out) :: status

        integer :: mode
        complex(dp) :: mode_rr, mode_rtheta, mode_rphi
        complex(dp) :: mode_thetatheta, mode_thetaphi, mode_phiphi

        radial_second = cmplx(0.0_dp, 0.0_dp, dp)
        radial_theta = radial_second
        radial_phi = radial_second
        theta_second = radial_second
        theta_phi = radial_second
        phi_second = radial_second
        if (.not. validate_expansion_inputs( &
                registry, coefficients, radius, theta, phi, status)) return
        do mode = 1, size(coefficients)
            call evaluate_mode_hessian( &
                registry, mode, radius, theta, phi, mode_rr, mode_rtheta, &
                mode_rphi, mode_thetatheta, mode_thetaphi, mode_phiphi, status)
            if (status%code /= FORTSPARSE_OK) return
            radial_second = radial_second + coefficients(mode)*mode_rr
            radial_theta = radial_theta + coefficients(mode)*mode_rtheta
            radial_phi = radial_phi + coefficients(mode)*mode_rphi
            theta_second = theta_second + coefficients(mode)*mode_thetatheta
            theta_phi = theta_phi + coefficients(mode)*mode_thetaphi
            phi_second = phi_second + coefficients(mode)*mode_phiphi
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fourier_mode_expansion_hessian

    subroutine evaluate_fourier_mode_expansion_hvp( &
            registry, coefficients, radius, theta, phi, radius_dot, theta_dot, &
            phi_dot, radial_gradient_dot, theta_gradient_dot, &
            phi_gradient_dot, status)
        !! Apply the fixed-coefficient coordinate Hessian to one direction.
        !!
        !! This is a deliberately scoped Hessian-vector product.  Coefficients
        !! remain fixed; coefficient/coordinate first-order tangents are already
        !! provided by `evaluate_fourier_mode_expansion_jvp`.
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: radius, theta, phi
        real(dp), intent(in) :: radius_dot, theta_dot, phi_dot
        complex(dp), intent(out) :: radial_gradient_dot
        complex(dp), intent(out) :: theta_gradient_dot, phi_gradient_dot
        type(fortsparse_status_t), intent(out) :: status

        complex(dp) :: radial_second, radial_theta, radial_phi
        complex(dp) :: theta_second, theta_phi, phi_second

        radial_gradient_dot = cmplx(0.0_dp, 0.0_dp, dp)
        theta_gradient_dot = radial_gradient_dot
        phi_gradient_dot = radial_gradient_dot
        if (.not. validate_expansion_inputs( &
                registry, coefficients, radius, theta, phi, status)) return
        if (.not. ieee_is_finite(radius_dot) .or. &
                .not. ieee_is_finite(theta_dot) .or. &
                .not. ieee_is_finite(phi_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier mode expansion HVP has invalid coordinate direction")
            return
        end if
        call evaluate_fourier_mode_expansion_hessian( &
            registry, coefficients, radius, theta, phi, radial_second, &
            radial_theta, radial_phi, theta_second, theta_phi, phi_second, status)
        if (status%code /= FORTSPARSE_OK) return
        radial_gradient_dot = radial_second*radius_dot + &
            radial_theta*theta_dot + radial_phi*phi_dot
        theta_gradient_dot = radial_theta*radius_dot + &
            theta_second*theta_dot + theta_phi*phi_dot
        phi_gradient_dot = radial_phi*radius_dot + theta_phi*theta_dot + &
            phi_second*phi_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fourier_mode_expansion_hvp

    subroutine evaluate_mode_hessian( &
            registry, mode, radius, theta, phi, radial_second, radial_theta, &
            radial_phi, theta_second, theta_phi, phi_second, status)
        type(fourier_mode_registry_t), intent(in) :: registry
        integer, intent(in) :: mode
        real(dp), intent(in) :: radius, theta, phi
        complex(dp), intent(out) :: radial_second, radial_theta, radial_phi
        complex(dp), intent(out) :: theta_second, theta_phi, phi_second
        type(fortsparse_status_t), intent(out) :: status

        complex(dp) :: mode_value, mode_radial, mode_theta, mode_phi
        complex(dp) :: phase_factor
        real(dp) :: phase, radial_factor_second
        real(dp) :: poloidal, toroidal

        call evaluate_fourier_mode( &
            registry, mode, radius, theta, phi, mode_value, mode_radial, &
            mode_theta, mode_phi, status)
        if (status%code /= FORTSPARSE_OK) return
        poloidal = real(registry%poloidal_modes(mode), dp)
        toroidal = real(registry%toroidal_modes(mode)*registry%field_periods, dp)
        phase = poloidal*(theta + registry%poloidal_phase) + &
            toroidal*(phi + registry%toroidal_phase)
        phase_factor = exp(cmplx(0.0_dp, phase, dp))
        radial_factor_second = radial_power_second( &
            radius, registry%radial_powers(mode))
        radial_second = registry%normalization(mode)*radial_factor_second* &
            phase_factor
        radial_theta = cmplx(0.0_dp, poloidal, dp)*mode_radial
        radial_phi = cmplx(0.0_dp, toroidal, dp)*mode_radial
        theta_second = cmplx(0.0_dp, poloidal, dp)*mode_theta
        theta_phi = cmplx(0.0_dp, toroidal, dp)*mode_theta
        phi_second = cmplx(0.0_dp, toroidal, dp)*mode_phi
    end subroutine evaluate_mode_hessian

    real(dp) function radial_power_second(radius, power) result(second)
        real(dp), intent(in) :: radius
        integer, intent(in) :: power

        if (power < 2) then
            second = 0.0_dp
        else
            second = real(power*(power - 1), dp)*radius**(power - 2)
        end if
    end function radial_power_second

    logical function validate_expansion_inputs( &
            registry, coefficients, radius, theta, phi, status)
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: radius, theta, phi
        type(fortsparse_status_t), intent(out) :: status

        validate_expansion_inputs = .false.
        if (.not. validate_fourier_mode_registry(registry, status)) return
        if (size(coefficients) /= size(registry%poloidal_modes) .or. &
                .not. finite_complex(coefficients) .or. &
                .not. ieee_is_finite(radius) .or. radius < 0.0_dp .or. &
                .not. ieee_is_finite(theta) .or. .not. ieee_is_finite(phi)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier mode expansion has invalid coefficients/coordinates")
            return
        end if
        validate_expansion_inputs = .true.
    end function validate_expansion_inputs

    logical function finite_complex(values) result(finite)
        complex(dp), intent(in) :: values(:)

        finite = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex

end module fortfem_fourier_mode_expansion
