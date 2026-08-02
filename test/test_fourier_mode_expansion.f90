program test_fourier_mode_expansion
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_fourier_mode_expansion, &
        evaluate_fourier_mode_expansion_jvp, &
        evaluate_fourier_mode_expansion_vjp
    use fortfem_fourier_mode_registry, only: &
        fourier_mode_registry_t, initialize_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: mode_count = 5
    integer, parameter :: poloidal(mode_count) = [-2, -1, 0, 1, 2]
    integer, parameter :: toroidal(mode_count) = [0, 1, 0, -1, 0]
    integer, parameter :: radial_power(mode_count) = [2, 1, 0, 1, 2]
    real(dp), parameter :: normalization(mode_count) = [ &
        0.8_dp, 1.1_dp, 0.7_dp, 1.3_dp, 0.9_dp]
    complex(dp), parameter :: coefficients(mode_count) = [ &
        cmplx(0.2_dp, -0.3_dp, dp), cmplx(-0.4_dp, 0.1_dp, dp), &
        cmplx(0.6_dp, 0.0_dp, dp), cmplx(0.3_dp, 0.5_dp, dp), &
        cmplx(-0.1_dp, 0.2_dp, dp)]
    complex(dp), parameter :: coefficients_dot(mode_count) = [ &
        cmplx(-0.05_dp, 0.07_dp, dp), cmplx(0.02_dp, -0.04_dp, dp), &
        cmplx(0.06_dp, 0.01_dp, dp), cmplx(-0.03_dp, 0.02_dp, dp), &
        cmplx(0.04_dp, -0.06_dp, dp)]
    complex(dp), parameter :: value_bar = cmplx(0.17_dp, -0.23_dp, dp)
    real(dp), parameter :: radius = 0.43_dp, theta = 0.71_dp, phi = -0.32_dp
    real(dp), parameter :: radius_dot = -0.08_dp, theta_dot = 0.12_dp
    real(dp), parameter :: phi_dot = 0.05_dp
    real(dp), parameter :: epsilon = 1.0e-7_dp
    real(dp), parameter :: tolerance = 2.0e-8_dp

    type(fourier_mode_registry_t) :: registry
    type(fortsparse_status_t) :: status
    complex(dp) :: value, radial_derivative, theta_derivative, phi_derivative
    complex(dp) :: value_dot, value_plus, value_minus
    complex(dp) :: expected, expected_radial, expected_theta, expected_phi
    complex(dp) :: expected_dot, coefficients_bar(mode_count)
    complex(dp) :: bad_coefficients(mode_count - 1), nonfinite_coefficients(mode_count)
    real(dp) :: radius_bar, theta_bar, phi_bar
    real(dp) :: lhs, rhs

    call initialize_fourier_mode_registry( &
        registry, poloidal, toroidal, 2, 0.13_dp, -0.21_dp, .false., &
        radial_powers=radial_power, normalization=normalization, status=status)
    call check_condition(status%code == 0, &
        "Fourier mode expansion registry initializes")

    call independent_expansion_oracle( &
        radius, theta, phi, expected, expected_radial, expected_theta, expected_phi)
    call evaluate_fourier_mode_expansion( &
        registry, coefficients, radius, theta, phi, value, radial_derivative, &
        theta_derivative, phi_derivative, status)
    call check_condition(status%code == 0 .and. &
        abs(value - expected) < tolerance .and. &
        abs(radial_derivative - expected_radial) < tolerance .and. &
        abs(theta_derivative - expected_theta) < tolerance .and. &
        abs(phi_derivative - expected_phi) < tolerance, &
        "Fourier mode expansion matches an independent phase/radial oracle")

    call independent_expansion_jvp_oracle( &
        radius, theta, phi, expected_dot)
    call evaluate_fourier_mode_expansion_jvp( &
        registry, coefficients, radius, theta, phi, coefficients_dot, &
        radius_dot, theta_dot, phi_dot, value_dot, status)
    call evaluate_fourier_mode_expansion( &
        registry, coefficients + epsilon*coefficients_dot, &
        radius + epsilon*radius_dot, theta + epsilon*theta_dot, &
        phi + epsilon*phi_dot, value_plus, radial_derivative, theta_derivative, &
        phi_derivative, status)
    call evaluate_fourier_mode_expansion( &
        registry, coefficients - epsilon*coefficients_dot, &
        radius - epsilon*radius_dot, theta - epsilon*theta_dot, &
        phi - epsilon*phi_dot, value_minus, radial_derivative, theta_derivative, &
        phi_derivative, status)
    call check_condition(status%code == 0 .and. &
        abs(value_dot - expected_dot) < tolerance .and. &
        abs((value_plus - value_minus)/(2.0_dp*epsilon) - value_dot) < tolerance, &
        "Fourier mode expansion JVP matches product and central differences")

    call evaluate_fourier_mode_expansion_vjp( &
        registry, coefficients, radius, theta, phi, value_bar, coefficients_bar, &
        radius_bar, theta_bar, phi_bar, status)
    lhs = real(conjg(value_bar)*expected_dot, dp)
    rhs = real(sum(conjg(coefficients_bar)*coefficients_dot), dp) + &
        radius_bar*radius_dot + theta_bar*theta_dot + phi_bar*phi_dot
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < tolerance, &
        "Fourier mode expansion VJP satisfies the real-part adjoint oracle")

    bad_coefficients = coefficients(:mode_count - 1)
    call evaluate_fourier_mode_expansion( &
        registry, bad_coefficients, radius, theta, phi, value, radial_derivative, &
        theta_derivative, phi_derivative, status)
    call check_condition(status%code /= 0, &
        "Fourier mode expansion rejects an incompatible coefficient shape")

    nonfinite_coefficients = coefficients
    nonfinite_coefficients(2) = cmplx( &
        ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    call evaluate_fourier_mode_expansion( &
        registry, nonfinite_coefficients, radius, theta, phi, value, &
        radial_derivative, theta_derivative, phi_derivative, status)
    call check_condition(status%code /= 0, &
        "Fourier mode expansion rejects non-finite coefficients")

    call evaluate_fourier_mode_expansion( &
        registry, coefficients, -0.1_dp, theta, phi, value, radial_derivative, &
        theta_derivative, phi_derivative, status)
    call check_condition(status%code /= 0, &
        "Fourier mode expansion rejects a negative radius")

    call check_summary("Fourier mode/radial expansion")

contains

    subroutine independent_expansion_oracle( &
            radius_value, theta_value, phi_value, value_oracle, radial_oracle, &
            theta_oracle, phi_oracle)
        real(dp), intent(in) :: radius_value, theta_value, phi_value
        complex(dp), intent(out) :: value_oracle, radial_oracle
        complex(dp), intent(out) :: theta_oracle, phi_oracle

        integer :: mode
        real(dp) :: radial_factor, radial_factor_derivative, phase
        complex(dp) :: phase_factor, basis, basis_radial

        value_oracle = cmplx(0.0_dp, 0.0_dp, dp)
        radial_oracle = value_oracle
        theta_oracle = value_oracle
        phi_oracle = value_oracle
        do mode = 1, mode_count
            radial_factor = radius_value**radial_power(mode)
            if (radial_power(mode) == 0) then
                radial_factor_derivative = 0.0_dp
            else
                radial_factor_derivative = real(radial_power(mode), dp)* &
                    radius_value**(radial_power(mode) - 1)
            end if
            phase = real(poloidal(mode), dp)*(theta_value + 0.13_dp) + &
                real(toroidal(mode), dp)*2.0_dp*(phi_value - 0.21_dp)
            phase_factor = cmplx(cos(phase), sin(phase), dp)
            basis = normalization(mode)*radial_factor*phase_factor
            basis_radial = normalization(mode)*radial_factor_derivative*phase_factor
            value_oracle = value_oracle + coefficients(mode)*basis
            radial_oracle = radial_oracle + coefficients(mode)*basis_radial
            theta_oracle = theta_oracle + coefficients(mode)* &
                cmplx(0.0_dp, real(poloidal(mode), dp), dp)*basis
            phi_oracle = phi_oracle + coefficients(mode)* &
                cmplx(0.0_dp, real(toroidal(mode)*2, dp), dp)*basis
        end do
    end subroutine independent_expansion_oracle

    subroutine independent_expansion_jvp_oracle( &
            radius_value, theta_value, phi_value, value_oracle)
        real(dp), intent(in) :: radius_value, theta_value, phi_value
        complex(dp), intent(out) :: value_oracle

        integer :: mode
        real(dp) :: phase, radial_factor, radial_factor_derivative
        complex(dp) :: phase_factor, basis, basis_radial
        complex(dp) :: basis_theta, basis_phi

        value_oracle = cmplx(0.0_dp, 0.0_dp, dp)
        do mode = 1, mode_count
            radial_factor = radius_value**radial_power(mode)
            if (radial_power(mode) == 0) then
                radial_factor_derivative = 0.0_dp
            else
                radial_factor_derivative = real(radial_power(mode), dp)* &
                    radius_value**(radial_power(mode) - 1)
            end if
            phase = real(poloidal(mode), dp)*(theta_value + 0.13_dp) + &
                real(toroidal(mode), dp)*2.0_dp*(phi_value - 0.21_dp)
            phase_factor = cmplx(cos(phase), sin(phase), dp)
            basis = normalization(mode)*radial_factor*phase_factor
            basis_radial = normalization(mode)*radial_factor_derivative*phase_factor
            basis_theta = cmplx(0.0_dp, real(poloidal(mode), dp), dp)*basis
            basis_phi = cmplx(0.0_dp, real(toroidal(mode)*2, dp), dp)*basis
            value_oracle = value_oracle + coefficients_dot(mode)*basis + &
                coefficients(mode)*(basis_radial*radius_dot + &
                basis_theta*theta_dot + basis_phi*phi_dot)
        end do
    end subroutine independent_expansion_jvp_oracle

end program test_fourier_mode_expansion
