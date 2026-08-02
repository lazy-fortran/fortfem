program test_toroidal_spectral_neumann_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_toroidal_spectral_trace, solve_toroidal_spectral_neumann, &
        solve_toroidal_spectral_neumann_jvp, solve_toroidal_spectral_neumann_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: mode_count = 3
    real(dp), parameter :: scale = 1.7_dp, eta = 1.1_dp, step = 1.0e-6_dp
    integer :: degree_indices(mode_count), orders(mode_count), status
    complex(dp) :: potential(mode_count), normal_coefficients(mode_count)
    complex(dp) :: potential_recovered(mode_count), potential_dot(mode_count)
    complex(dp) :: normal_coefficients_dot(mode_count)
    complex(dp) :: potential_plus(mode_count), potential_minus(mode_count)
    complex(dp) :: normal_coefficients_bar(mode_count), potential_bar(mode_count)
    real(dp) :: scale_dot, eta_dot, scale_bar, eta_bar
    real(dp) :: derivative_error, adjoint_error, lhs, rhs
    complex(dp) :: value, normal
    logical :: all_passed
    integer :: mode

    all_passed = .true.
    degree_indices = [0, 2, 3]
    orders = [0, 1, 2]
    potential = [ &
        cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.2_dp, 0.6_dp, dp)]
    normal_coefficients_dot = [ &
        cmplx(0.03_dp, 0.01_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp), &
        cmplx(0.01_dp, -0.03_dp, dp)]
    scale_dot = 0.017_dp
    eta_dot = -0.013_dp

    do mode = 1, mode_count
        call evaluate_toroidal_spectral_trace( &
            degree_indices(mode:mode), orders(mode:mode), potential(mode:mode), &
            scale, eta, 0.0_dp, 0.0_dp, .false., value, normal, status)
        if (status == 0) normal_coefficients(mode) = normal
    end do
    call solve_toroidal_spectral_neumann( &
        degree_indices, orders, normal_coefficients, scale, eta, .false., .true., &
        1.0e-12_dp, potential_recovered, status)
    call record_condition(status == 0 .and. maxval(abs(potential_recovered - potential)) < &
        2.0e-13_dp, "toroidal modal Neumann solve recovers supplied coefficients")
    call solve_toroidal_spectral_neumann( &
        degree_indices, orders, normal_coefficients, scale, eta, .false., .false., &
        1.0e-12_dp, potential_recovered, status)
    call record_condition(status /= 0, "toroidal modal Neumann solve exposes zero-mode policy")

    call solve_toroidal_spectral_neumann_jvp( &
        degree_indices, orders, normal_coefficients, scale, eta, .false., .true., &
        1.0e-12_dp, normal_coefficients_dot, scale_dot, eta_dot, potential_dot, status)
    call solve_toroidal_spectral_neumann( &
        degree_indices, orders, normal_coefficients + step*normal_coefficients_dot, &
        scale + step*scale_dot, eta + step*eta_dot, .false., .true., 1.0e-12_dp, &
        potential_plus, status)
    call solve_toroidal_spectral_neumann( &
        degree_indices, orders, normal_coefficients - step*normal_coefficients_dot, &
        scale - step*scale_dot, eta - step*eta_dot, .false., .true., 1.0e-12_dp, &
        potential_minus, status)
    derivative_error = maxval(abs(potential_dot - &
        (potential_plus - potential_minus)/(2.0_dp*step)))
    call record_condition(status == 0 .and. derivative_error < 2.0e-8_dp, &
        "toroidal modal Neumann JVP matches central reassembly")

    potential_bar = [ &
        cmplx(0.8_dp, -0.3_dp, dp), cmplx(-0.2_dp, 0.6_dp, dp), &
        cmplx(0.1_dp, 0.4_dp, dp)]
    call solve_toroidal_spectral_neumann_vjp( &
        degree_indices, orders, normal_coefficients, scale, eta, .false., .true., &
        1.0e-12_dp, potential_bar, normal_coefficients_bar, scale_bar, eta_bar, status)
    lhs = real(sum(conjg(potential_bar)*potential_dot), dp)
    rhs = real(sum(conjg(normal_coefficients_bar)*normal_coefficients_dot), dp) + &
        scale_bar*scale_dot + eta_bar*eta_dot
    adjoint_error = abs(lhs - rhs)
    call record_condition(status == 0 .and. adjoint_error < &
        5.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "toroidal modal Neumann VJP satisfies the real complex adjoint")

    call check_summary("toroidal modal Neumann")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_toroidal_spectral_neumann_ad
