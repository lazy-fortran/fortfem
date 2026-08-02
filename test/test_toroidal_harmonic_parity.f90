program test_toroidal_harmonic_parity
    !! Independent phase/parity oracle for toroidal P/Q traces.
    !!
    !! The spectral trace implementation is exercised against the separated
    !! expression written here directly.  This keeps the Fourier convention
    !! and the outward normal sign visible at the boundary between the
    !! special-function and FEM/BEM layers.
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_toroidal_spectral_trace, toroidal_p, &
        toroidal_q, toroidal_p_derivative, toroidal_q_derivative
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: degree_index = 3, order = 2
    real(dp), parameter :: scale = 1.4_dp, eta = 1.15_dp
    real(dp), parameter :: theta = 0.41_dp, phi = -0.73_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: radial_step = 2.0e-4_dp
    real(dp), parameter :: tolerance = 2.0e-10_dp
    integer :: status
    complex(dp) :: value, normal, value_periodic, normal_periodic
    complex(dp) :: value_reflected, normal_reflected
    complex(dp) :: expected_value, expected_normal
    complex(dp) :: coefficient, coefficient_sum, normal_sum
    real(dp) :: ode_residual_p, ode_residual_q
    logical :: all_passed

    all_passed = .true.
    coefficient = cmplx(0.7_dp, -0.2_dp, dp)

    call check_trace(.false., "P")
    call check_trace(.true., "Q")

    ! A real radial branch with a real coefficient is conjugate-reflected by
    ! (theta,phi)->(-theta,-phi), while both angular periods are exact.
    call evaluate_toroidal_spectral_trace( &
        [degree_index], [order], [cmplx(1.0_dp, 0.0_dp, dp)], scale, eta, &
        theta, phi, .false., value, normal, status)
    call evaluate_toroidal_spectral_trace( &
        [degree_index], [order], [cmplx(1.0_dp, 0.0_dp, dp)], scale, eta, &
        theta + 2.0_dp*pi, phi + 2.0_dp*pi, .false., value_periodic, &
        normal_periodic, status)
    call record(status == 0 .and. max(abs(value_periodic - value), &
        abs(normal_periodic - normal)) < tolerance, &
        "toroidal P trace is periodic in both angles")
    call evaluate_toroidal_spectral_trace( &
        [degree_index], [order], [cmplx(1.0_dp, 0.0_dp, dp)], scale, eta, &
        -theta, -phi, .false., value_reflected, normal_reflected, status)
    call record(status == 0 .and. max(abs(value_reflected - conjg(value)), &
        abs(normal_reflected - conjg(normal))) < tolerance, &
        "real P trace has the expected reflection parity")

    call direct_trace(.false., coefficient, expected_value, expected_normal)
    call evaluate_toroidal_spectral_trace( &
        [degree_index], [order], [coefficient], scale, eta, theta, phi, .false., &
        value, normal, status)
    call record(status == 0 .and. max(abs(value - expected_value), &
        abs(normal - expected_normal)) < tolerance, &
        "P trace matches the direct separated normalization and normal sign")

    call direct_trace(.true., coefficient, expected_value, expected_normal)
    call evaluate_toroidal_spectral_trace( &
        [degree_index], [order], [coefficient], scale, eta, theta, phi, .true., &
        value, normal, status)
    call record(status == 0 .and. max(abs(value - expected_value), &
        abs(normal - expected_normal)) < tolerance, &
        "Q trace matches the direct separated normalization and normal sign")

    coefficient_sum = cmplx(0.0_dp, 0.0_dp, dp)
    normal_sum = cmplx(0.0_dp, 0.0_dp, dp)
    call direct_trace(.false., cmplx(1.0_dp, 0.0_dp, dp), expected_value, &
        expected_normal)
    coefficient_sum = coefficient_sum + expected_value
    normal_sum = normal_sum + expected_normal
    call direct_trace(.false., cmplx(-0.3_dp, 0.4_dp, dp), expected_value, &
        expected_normal)
    coefficient_sum = coefficient_sum + expected_value
    normal_sum = normal_sum + expected_normal
    call evaluate_toroidal_spectral_trace( &
        [degree_index, degree_index], [order, order], &
        [cmplx(1.0_dp, 0.0_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp)], scale, eta, &
        theta, phi, .false., value, normal, status)
    call record(status == 0 .and. abs(value - coefficient_sum) < tolerance .and. &
        abs(normal - normal_sum) < tolerance, &
        "modal coefficients contract the normalized P trace linearly")

    call record(associated_ode_residual(.false., ode_residual_p) .and. &
        ode_residual_p < 2.0e-5_dp, &
        "P radial branch satisfies the associated Legendre ODE")
    call record(associated_ode_residual(.true., ode_residual_q) .and. &
        ode_residual_q < 2.0e-5_dp, &
        "Q radial branch satisfies the associated Legendre ODE")

    call check_summary("Toroidal harmonic parity")
    if (.not. all_passed) error stop 1

contains

    subroutine check_trace(use_second_kind, label)
        logical, intent(in) :: use_second_kind
        character(*), intent(in) :: label
        real(dp) :: radial, radial_first, denominator, square_root

        if (use_second_kind) then
            radial = toroidal_q(degree_index, order, cosh(eta))
            radial_first = toroidal_q_derivative( &
                degree_index, order, cosh(eta))
        else
            radial = toroidal_p(degree_index, order, cosh(eta))
            radial_first = toroidal_p_derivative( &
                degree_index, order, cosh(eta))
        end if
        denominator = cosh(eta) - cos(theta)
        square_root = sqrt(denominator)
        expected_value = square_root*radial*exp(cmplx(0.0_dp, &
            real(degree_index, dp)*theta + real(order, dp)*phi, dp))
        expected_normal = -denominator/scale*sinh(eta)* &
            (radial/(2.0_dp*square_root) + square_root*radial_first)* &
            exp(cmplx(0.0_dp, real(degree_index, dp)*theta + &
            real(order, dp)*phi, dp))
        call evaluate_toroidal_spectral_trace( &
            [degree_index], [order], [cmplx(1.0_dp, 0.0_dp, dp)], scale, eta, &
            theta, phi, use_second_kind, value, normal, status)
        call record(status == 0 .and. max(abs(value - expected_value), &
            abs(normal - expected_normal)) < tolerance, &
            label//" trace follows the direct separated expression")
    end subroutine check_trace

    subroutine direct_trace(use_second_kind, factor, value_out, normal_out)
        logical, intent(in) :: use_second_kind
        complex(dp), intent(in) :: factor
        complex(dp), intent(out) :: value_out, normal_out
        real(dp) :: radial, radial_first, denominator, square_root
        complex(dp) :: phase

        if (use_second_kind) then
            radial = toroidal_q(degree_index, order, cosh(eta))
            radial_first = toroidal_q_derivative( &
                degree_index, order, cosh(eta))
        else
            radial = toroidal_p(degree_index, order, cosh(eta))
            radial_first = toroidal_p_derivative( &
                degree_index, order, cosh(eta))
        end if
        denominator = cosh(eta) - cos(theta)
        square_root = sqrt(denominator)
        phase = exp(cmplx(0.0_dp, real(degree_index, dp)*theta + &
            real(order, dp)*phi, dp))
        value_out = factor*square_root*radial*phase
        normal_out = factor*(-denominator/scale)*sinh(eta)* &
            (radial/(2.0_dp*square_root) + square_root*radial_first)*phase
    end subroutine direct_trace

    logical function associated_ode_residual(use_second_kind, residual) result(ok)
        logical, intent(in) :: use_second_kind
        real(dp), intent(out) :: residual
        real(dp) :: x, radial_minus, radial_zero, radial_plus
        real(dp) :: first_derivative, second_derivative, degree

        x = cosh(eta)
        if (use_second_kind) then
            radial_minus = toroidal_q(degree_index, order, x - radial_step)
            radial_zero = toroidal_q(degree_index, order, x)
            radial_plus = toroidal_q(degree_index, order, x + radial_step)
        else
            radial_minus = toroidal_p(degree_index, order, x - radial_step)
            radial_zero = toroidal_p(degree_index, order, x)
            radial_plus = toroidal_p(degree_index, order, x + radial_step)
        end if
        first_derivative = (radial_plus - radial_minus)/(2.0_dp*radial_step)
        second_derivative = (radial_plus - 2.0_dp*radial_zero + radial_minus)/ &
            radial_step**2
        degree = real(degree_index, dp) - 0.5_dp
        residual = abs((x*x - 1.0_dp)*second_derivative + 2.0_dp*x* &
            first_derivative - (degree*(degree + 1.0_dp) + &
            real(order*order, dp)/(x*x - 1.0_dp))*radial_zero)
        ok = residual == residual
    end function associated_ode_residual

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_toroidal_harmonic_parity
