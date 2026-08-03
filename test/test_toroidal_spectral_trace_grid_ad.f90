program test_toroidal_spectral_trace_grid_ad
    use check, only: check_condition, check_summary
    use fortfem_fourier, only: &
        evaluate_toroidal_spectral_trace, evaluate_toroidal_spectral_trace_grid, &
        evaluate_toroidal_spectral_trace_grid_jvp, &
        evaluate_toroidal_spectral_trace_grid_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: point_count = 3
    real(dp), parameter :: scale = 1.7_dp, step = 1.0e-5_dp
    integer :: degree_indices(3), orders(3), status, point
    complex(dp) :: coefficients(3), coefficients_dot(3)
    complex(dp), allocatable :: coefficients_bar(:)
    real(dp) :: eta(point_count), theta(point_count), phi(point_count)
    real(dp) :: eta_dot(point_count), theta_dot(point_count), phi_dot(point_count)
    real(dp) :: eta_bar(point_count), theta_bar(point_count), phi_bar(point_count)
    real(dp) :: scale_dot, scale_bar
    complex(dp) :: values(point_count), normal_derivatives(point_count)
    complex(dp) :: values_dot(point_count), normal_derivatives_dot(point_count)
    complex(dp) :: values_plus(point_count), values_minus(point_count)
    complex(dp) :: normals_plus(point_count), normals_minus(point_count)
    complex(dp) :: values_bar(point_count), normal_derivatives_bar(point_count)
    complex(dp) :: value_single, normal_single
    real(dp) :: derivative_error, adjoint_error, lhs, rhs
    logical :: all_passed

    all_passed = .true.
    degree_indices = [1, 2, 3]
    orders = [0, 1, 2]
    coefficients = [ &
        cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.2_dp, 0.6_dp, dp)]
    coefficients_dot = [ &
        cmplx(0.03_dp, 0.01_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp), &
        cmplx(0.01_dp, -0.03_dp, dp)]
    eta = [0.8_dp, 1.1_dp, 1.4_dp]
    theta = [0.2_dp, 1.1_dp, 2.4_dp]
    phi = [-0.6_dp, 0.3_dp, 1.7_dp]
    eta_dot = [-0.013_dp, 0.017_dp, -0.009_dp]
    theta_dot = [0.021_dp, -0.015_dp, 0.011_dp]
    phi_dot = [-0.019_dp, 0.014_dp, 0.008_dp]
    scale_dot = 0.017_dp

    call evaluate_toroidal_spectral_trace_grid( &
        degree_indices, orders, coefficients, scale, eta, theta, phi, .false., &
        values, normal_derivatives, status)
    call record_condition(status == 0, &
        "toroidal spectral trace grid accepts periodic samples")
    do point = 1, point_count
        call evaluate_toroidal_spectral_trace( &
            degree_indices, orders, coefficients, scale, eta(point), theta(point), &
            phi(point), .false., value_single, normal_single, status)
        call record_condition(status == 0 .and. &
            abs(values(point) - value_single) < 2.0e-13_dp .and. &
            abs(normal_derivatives(point) - normal_single) < 2.0e-13_dp, &
            "grid trace agrees with the scalar modal trace")
    end do

    call evaluate_toroidal_spectral_trace_grid_jvp( &
        degree_indices, orders, coefficients, scale, eta, theta, phi, .false., &
        coefficients_dot, scale_dot, eta_dot, theta_dot, phi_dot, values_dot, &
        normal_derivatives_dot, status)
    call evaluate_toroidal_spectral_trace_grid( &
        degree_indices, orders, coefficients + step*coefficients_dot, &
        scale + step*scale_dot, eta + step*eta_dot, theta + step*theta_dot, &
        phi + step*phi_dot, .false., values_plus, normals_plus, status)
    call evaluate_toroidal_spectral_trace_grid( &
        degree_indices, orders, coefficients - step*coefficients_dot, &
        scale - step*scale_dot, eta - step*eta_dot, theta - step*theta_dot, &
        phi - step*phi_dot, .false., values_minus, normals_minus, status)
    derivative_error = max( &
        maxval(abs(values_dot - (values_plus - values_minus)/(2.0_dp*step))), &
        maxval(abs(normal_derivatives_dot - &
        (normals_plus - normals_minus)/(2.0_dp*step))))
    call record_condition(status == 0 .and. derivative_error < 2.0e-7_dp, &
        "toroidal spectral trace grid JVP matches central reassembly")

    values_bar = [ &
        cmplx(0.8_dp, -0.3_dp, dp), cmplx(-0.2_dp, 0.6_dp, dp), &
        cmplx(0.4_dp, 0.1_dp, dp)]
    normal_derivatives_bar = [ &
        cmplx(-0.3_dp, 0.2_dp, dp), cmplx(0.5_dp, -0.1_dp, dp), &
        cmplx(-0.4_dp, 0.3_dp, dp)]
    call evaluate_toroidal_spectral_trace_grid_vjp( &
        degree_indices, orders, coefficients, scale, eta, theta, phi, .false., &
        values_bar, normal_derivatives_bar, coefficients_bar, scale_bar, eta_bar, &
        theta_bar, phi_bar, status)
    lhs = real(sum(conjg(values_bar)*values_dot) + &
        sum(conjg(normal_derivatives_bar)*normal_derivatives_dot), dp)
    rhs = real(sum(conjg(coefficients_bar)*coefficients_dot), dp) + &
        scale_bar*scale_dot + dot_product(eta_bar, eta_dot) + &
        dot_product(theta_bar, theta_dot) + dot_product(phi_bar, phi_dot)
    adjoint_error = abs(lhs - rhs)
    call record_condition(status == 0 .and. adjoint_error < &
        5.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "toroidal spectral trace grid VJP satisfies the real complex adjoint")

    call check_summary("toroidal spectral trace grid")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_toroidal_spectral_trace_grid_ad
