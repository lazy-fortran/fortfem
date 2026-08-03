program test_fourier_mode_expansion_hessian
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use check, only: check_condition, check_summary
    use fortfem_fourier, only: &
        evaluate_fourier_mode_expansion, &
        evaluate_fourier_mode_expansion_hessian, &
        evaluate_fourier_mode_expansion_hvp, &
        fourier_mode_registry_t, initialize_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: mode_count = 4
    integer, parameter :: poloidal(mode_count) = [-1, 0, 1, 2]
    integer, parameter :: toroidal(mode_count) = [1, 0, -1, 0]
    integer, parameter :: radial_power(mode_count) = [1, 0, 1, 2]
    real(dp), parameter :: normalization(mode_count) = [ &
        1.1_dp, 0.8_dp, 0.9_dp, 1.2_dp]
    complex(dp), parameter :: coefficients(mode_count) = [ &
        cmplx(0.2_dp, -0.3_dp, dp), cmplx(-0.4_dp, 0.1_dp, dp), &
        cmplx(0.6_dp, 0.2_dp, dp), cmplx(-0.1_dp, 0.5_dp, dp)]
    real(dp), parameter :: radius = 0.47_dp, theta = 0.63_dp, phi = -0.29_dp
    real(dp), parameter :: radius_dot = -0.09_dp, theta_dot = 0.14_dp
    real(dp), parameter :: phi_dot = 0.06_dp
    real(dp), parameter :: epsilon = 1.0e-7_dp
    real(dp), parameter :: tolerance = 3.0e-8_dp

    type(fourier_mode_registry_t) :: registry
    type(fortsparse_status_t) :: status
    complex(dp) :: radial_second, radial_theta, radial_phi
    complex(dp) :: theta_second, theta_phi, phi_second
    complex(dp) :: expected_rr, expected_rt, expected_rp
    complex(dp) :: expected_tt, expected_tp, expected_pp
    complex(dp) :: radial_gradient_dot, theta_gradient_dot, phi_gradient_dot
    complex(dp) :: value_plus, value_minus
    complex(dp) :: radial_plus, theta_plus, phi_plus
    complex(dp) :: radial_minus, theta_minus, phi_minus
    complex(dp) :: hvp_fd(3)

    call initialize_fourier_mode_registry( &
        registry, poloidal, toroidal, 2, 0.11_dp, -0.17_dp, .false., &
        radial_powers=radial_power, normalization=normalization, status=status)
    call check_condition(status%code == 0, &
        "Fourier Hessian registry initializes")

    call independent_hessian_oracle( &
        expected_rr, expected_rt, expected_rp, expected_tt, expected_tp, expected_pp)
    call evaluate_fourier_mode_expansion_hessian( &
        registry, coefficients, radius, theta, phi, radial_second, radial_theta, &
        radial_phi, theta_second, theta_phi, phi_second, status)
    call check_condition(status%code == 0 .and. &
        abs(radial_second - expected_rr) < tolerance .and. &
        abs(radial_theta - expected_rt) < tolerance .and. &
        abs(radial_phi - expected_rp) < tolerance .and. &
        abs(theta_second - expected_tt) < tolerance .and. &
        abs(theta_phi - expected_tp) < tolerance .and. &
        abs(phi_second - expected_pp) < tolerance, &
        "Fourier expansion Hessian matches independent nested-loop oracle")

    call evaluate_fourier_mode_expansion_hvp( &
        registry, coefficients, radius, theta, phi, radius_dot, theta_dot, &
        phi_dot, radial_gradient_dot, theta_gradient_dot, phi_gradient_dot, status)
    call evaluate_fourier_mode_expansion( &
        registry, coefficients, radius + epsilon*radius_dot, &
        theta + epsilon*theta_dot, phi + epsilon*phi_dot, value_plus, &
        radial_plus, theta_plus, phi_plus, status)
    call evaluate_fourier_mode_expansion( &
        registry, coefficients, radius - epsilon*radius_dot, &
        theta - epsilon*theta_dot, phi - epsilon*phi_dot, value_minus, &
        radial_minus, theta_minus, phi_minus, status)
    hvp_fd = [(radial_plus - radial_minus)/(2.0_dp*epsilon), &
        (theta_plus - theta_minus)/(2.0_dp*epsilon), &
        (phi_plus - phi_minus)/(2.0_dp*epsilon)]
    call check_condition(status%code == 0 .and. &
        abs(radial_gradient_dot - hvp_fd(1)) < tolerance .and. &
        abs(theta_gradient_dot - hvp_fd(2)) < tolerance .and. &
        abs(phi_gradient_dot - hvp_fd(3)) < tolerance, &
        "Fourier expansion HVP matches finite differences of the gradient")

    call evaluate_fourier_mode_expansion_hvp( &
        registry, coefficients, radius, theta, phi, radius_dot, theta_dot, &
        ieee_value(0.0_dp, ieee_quiet_nan), radial_gradient_dot, theta_gradient_dot, &
        phi_gradient_dot, status)
    call check_condition(status%code /= 0, &
        "Fourier expansion HVP rejects a non-finite direction")

    call check_summary("Fourier mode/radial Hessian")

contains

    subroutine independent_hessian_oracle( &
            rr, rt, rp, tt, tp, pp)
        complex(dp), intent(out) :: rr, rt, rp, tt, tp, pp

        integer :: mode
        real(dp) :: phase, factor, first, second, poloidal_value, toroidal_value
        complex(dp) :: phase_factor, basis, basis_r, basis_rr
        complex(dp) :: basis_t, basis_p

        rr = cmplx(0.0_dp, 0.0_dp, dp)
        rt = rr
        rp = rr
        tt = rr
        tp = rr
        pp = rr
        do mode = 1, mode_count
            poloidal_value = real(poloidal(mode), dp)
            toroidal_value = real(2*toroidal(mode), dp)
            phase = poloidal_value*(theta + 0.11_dp) + &
                toroidal_value*(phi - 0.17_dp)
            phase_factor = cmplx(cos(phase), sin(phase), dp)
            factor = radius**radial_power(mode)
            if (radial_power(mode) == 0) then
                first = 0.0_dp
                second = 0.0_dp
            else
                first = real(radial_power(mode), dp)* &
                    radius**(radial_power(mode) - 1)
                if (radial_power(mode) == 1) then
                    second = 0.0_dp
                else
                    second = real(radial_power(mode)*(radial_power(mode) - 1), dp)* &
                        radius**(radial_power(mode) - 2)
                end if
            end if
            basis = normalization(mode)*factor*phase_factor
            basis_r = normalization(mode)*first*phase_factor
            basis_rr = normalization(mode)*second*phase_factor
            basis_t = cmplx(0.0_dp, poloidal_value, dp)*basis
            basis_p = cmplx(0.0_dp, toroidal_value, dp)*basis
            rr = rr + coefficients(mode)*basis_rr
            rt = rt + coefficients(mode)*cmplx(0.0_dp, poloidal_value, dp)*basis_r
            rp = rp + coefficients(mode)*cmplx(0.0_dp, toroidal_value, dp)*basis_r
            tt = tt + coefficients(mode)*cmplx(0.0_dp, poloidal_value, dp)*basis_t
            tp = tp + coefficients(mode)*cmplx(0.0_dp, toroidal_value, dp)*basis_t
            pp = pp + coefficients(mode)*cmplx(0.0_dp, toroidal_value, dp)*basis_p
        end do
    end subroutine independent_hessian_oracle

end program test_fourier_mode_expansion_hessian
