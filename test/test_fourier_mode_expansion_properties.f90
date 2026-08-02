program test_fourier_mode_expansion_properties
    use, intrinsic :: iso_fortran_env, only: int32
    use check, only: check_condition, check_property, check_summary, &
        property_random_integer, property_random_unit, property_rng_initialize, &
        property_rng_t
    use fortfem_api, only: &
        evaluate_fourier_mode_expansion, &
        evaluate_fourier_mode_expansion_hessian, &
        evaluate_fourier_mode_expansion_hvp, &
        fourier_mode_registry_t, initialize_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    logical :: all_passed
    integer(int32) :: first_failed_seed, shrunk_seed

    call check_property( &
        "random Fourier/radial expansions preserve value and Hessian algebra", &
        424242_int32, 24, mode_case, all_passed, first_failed_seed, shrunk_seed)
    call check_condition(all_passed .and. first_failed_seed == 0_int32 .and. &
        shrunk_seed == 0_int32, &
        "random Fourier expansion property reports no failure seed")
    call check_summary("random Fourier mode expansion properties")
    if (.not. all_passed) error stop 1

contains

    logical function mode_case(case_seed)
        integer(int32), intent(in) :: case_seed
        integer, parameter :: mode_count = 5
        type(property_rng_t) :: rng
        type(fourier_mode_registry_t) :: registry
        type(fortsparse_status_t) :: status
        integer :: poloidal(mode_count), toroidal(mode_count)
        integer :: radial_power(mode_count), mode, field_periods
        real(dp) :: normalization(mode_count), poloidal_phase, toroidal_phase
        real(dp) :: radius, theta, phi, radius_dot, theta_dot, phi_dot
        complex(dp) :: coefficients(mode_count)
        complex(dp) :: value, radial_derivative, theta_derivative, phi_derivative
        complex(dp) :: radial_second, radial_theta, radial_phi
        complex(dp) :: theta_second, theta_phi, phi_second
        complex(dp) :: radial_gradient_dot, theta_gradient_dot, phi_gradient_dot
        complex(dp) :: expected_value, expected_gradient(3), expected_hessian(6)
        complex(dp) :: expected_hvp(3)
        real(dp), parameter :: tolerance = 3.0e-11_dp

        mode_case = .false.
        call property_rng_initialize(rng, case_seed)
        poloidal = [-2, -1, 0, 1, 2]
        toroidal = [0, 1, 0, -1, 0]
        radial_power(1) = 2 + 2*property_random_integer(rng, 0, 1)
        radial_power(2) = 1 + 2*property_random_integer(rng, 0, 1)
        radial_power(3) = 2*property_random_integer(rng, 0, 2)
        radial_power(4) = 1 + 2*property_random_integer(rng, 0, 1)
        radial_power(5) = 2 + 2*property_random_integer(rng, 0, 1)
        do mode = 1, mode_count
            normalization(mode) = 0.4_dp + property_random_unit(rng)
            coefficients(mode) = cmplx( &
                2.0_dp*property_random_unit(rng) - 1.0_dp, &
                2.0_dp*property_random_unit(rng) - 1.0_dp, dp)
        end do
        field_periods = property_random_integer(rng, 1, 3)
        poloidal_phase = property_random_unit(rng) - 0.5_dp
        toroidal_phase = property_random_unit(rng) - 0.5_dp
        radius = 0.1_dp + 0.8_dp*property_random_unit(rng)
        theta = 2.0_dp*acos(-1.0_dp)*property_random_unit(rng) - acos(-1.0_dp)
        phi = 2.0_dp*acos(-1.0_dp)*property_random_unit(rng) - acos(-1.0_dp)
        radius_dot = 0.6_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
        theta_dot = 0.6_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
        phi_dot = 0.6_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)

        call initialize_fourier_mode_registry( &
            registry, poloidal, toroidal, field_periods, poloidal_phase, &
            toroidal_phase, .false., radial_power, normalization, status)
        if (status%code /= 0) return
        call independent_oracle( &
            poloidal, toroidal, radial_power, normalization, coefficients, &
            field_periods, poloidal_phase, toroidal_phase, radius, theta, phi, &
            radius_dot, theta_dot, phi_dot, expected_value, expected_gradient, &
            expected_hessian, expected_hvp)

        call evaluate_fourier_mode_expansion( &
            registry, coefficients, radius, theta, phi, value, radial_derivative, &
            theta_derivative, phi_derivative, status)
        if (status%code /= 0 .or. abs(value - expected_value) > tolerance .or. &
                abs(radial_derivative - expected_gradient(1)) > tolerance .or. &
                abs(theta_derivative - expected_gradient(2)) > tolerance .or. &
                abs(phi_derivative - expected_gradient(3)) > tolerance) return

        call evaluate_fourier_mode_expansion_hessian( &
            registry, coefficients, radius, theta, phi, radial_second, &
            radial_theta, radial_phi, theta_second, theta_phi, phi_second, status)
        if (status%code /= 0 .or. &
                abs(radial_second - expected_hessian(1)) > tolerance .or. &
                abs(radial_theta - expected_hessian(2)) > tolerance .or. &
                abs(radial_phi - expected_hessian(3)) > tolerance .or. &
                abs(theta_second - expected_hessian(4)) > tolerance .or. &
                abs(theta_phi - expected_hessian(5)) > tolerance .or. &
                abs(phi_second - expected_hessian(6)) > tolerance) return

        call evaluate_fourier_mode_expansion_hvp( &
            registry, coefficients, radius, theta, phi, radius_dot, theta_dot, &
            phi_dot, radial_gradient_dot, theta_gradient_dot, phi_gradient_dot, &
            status)
        if (status%code /= 0 .or. &
                abs(radial_gradient_dot - expected_hvp(1)) > tolerance .or. &
                abs(theta_gradient_dot - expected_hvp(2)) > tolerance .or. &
                abs(phi_gradient_dot - expected_hvp(3)) > tolerance) return
        mode_case = .true.
    end function mode_case

    subroutine independent_oracle( &
            poloidal, toroidal, radial_power, normalization, coefficients, &
            field_periods, poloidal_phase, toroidal_phase, radius, theta, phi, &
            radius_dot, theta_dot, phi_dot, value, gradient, hessian, hvp)
        integer, intent(in) :: poloidal(:), toroidal(:), radial_power(:)
        real(dp), intent(in) :: normalization(:)
        complex(dp), intent(in) :: coefficients(:)
        integer, intent(in) :: field_periods
        real(dp), intent(in) :: poloidal_phase, toroidal_phase
        real(dp), intent(in) :: radius, theta, phi, radius_dot, theta_dot, phi_dot
        complex(dp), intent(out) :: value, gradient(:), hessian(:), hvp(:)

        integer :: mode, power
        real(dp) :: phase, factor, first_factor, second_factor
        real(dp) :: poloidal_value, toroidal_value
        complex(dp) :: phase_factor, basis, basis_r, basis_rr
        complex(dp) :: basis_t, basis_p

        value = cmplx(0.0_dp, 0.0_dp, dp)
        gradient = value
        hessian = value
        do mode = 1, size(coefficients)
            power = radial_power(mode)
            poloidal_value = real(poloidal(mode), dp)
            toroidal_value = real(toroidal(mode)*field_periods, dp)
            phase = poloidal_value*(theta + poloidal_phase) + &
                toroidal_value*(phi + toroidal_phase)
            phase_factor = cmplx(cos(phase), sin(phase), dp)
            factor = radius**power
            if (power == 0) then
                first_factor = 0.0_dp
                second_factor = 0.0_dp
            else
                first_factor = real(power, dp)*radius**(power - 1)
                if (power == 1) then
                    second_factor = 0.0_dp
                else
                    second_factor = real(power*(power - 1), dp)* &
                        radius**(power - 2)
                end if
            end if
            basis = normalization(mode)*factor*phase_factor
            basis_r = normalization(mode)*first_factor*phase_factor
            basis_rr = normalization(mode)*second_factor*phase_factor
            basis_t = cmplx(0.0_dp, poloidal_value, dp)*basis
            basis_p = cmplx(0.0_dp, toroidal_value, dp)*basis
            value = value + coefficients(mode)*basis
            gradient(1) = gradient(1) + coefficients(mode)*basis_r
            gradient(2) = gradient(2) + coefficients(mode)*basis_t
            gradient(3) = gradient(3) + coefficients(mode)*basis_p
            hessian(1) = hessian(1) + coefficients(mode)*basis_rr
            hessian(2) = hessian(2) + coefficients(mode)* &
                cmplx(0.0_dp, poloidal_value, dp)*basis_r
            hessian(3) = hessian(3) + coefficients(mode)* &
                cmplx(0.0_dp, toroidal_value, dp)*basis_r
            hessian(4) = hessian(4) + coefficients(mode)* &
                cmplx(0.0_dp, poloidal_value, dp)*basis_t
            hessian(5) = hessian(5) + coefficients(mode)* &
                cmplx(0.0_dp, toroidal_value, dp)*basis_t
            hessian(6) = hessian(6) + coefficients(mode)* &
                cmplx(0.0_dp, toroidal_value, dp)*basis_p
        end do
        hvp(1) = hessian(1)*radius_dot + hessian(2)*theta_dot + &
            hessian(3)*phi_dot
        hvp(2) = hessian(2)*radius_dot + hessian(4)*theta_dot + &
            hessian(5)*phi_dot
        hvp(3) = hessian(3)*radius_dot + hessian(5)*theta_dot + &
            hessian(6)*phi_dot
    end subroutine independent_oracle

end program test_fourier_mode_expansion_properties
