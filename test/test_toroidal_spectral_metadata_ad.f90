program test_toroidal_spectral_metadata_ad
    use check, only: check_condition, check_summary
    use fortfem_fourier, only: analyze_toroidal_spectral_modes, &
        analyze_toroidal_spectral_modes_jvp, analyze_toroidal_spectral_modes_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: mode_count = 4
    real(dp), parameter :: step = 1.0e-6_dp
    integer :: degree_indices(mode_count), orders(mode_count)
    integer :: degree_limit, order_limit, retained_count, omitted_count
    integer :: zero_mode_count, status
    complex(dp) :: coefficients(mode_count), coefficients_dot(mode_count)
    complex(dp) :: coefficients_bar(mode_count)
    real(dp) :: total_energy, omitted_energy, total_energy_dot, omitted_energy_dot
    real(dp) :: total_plus, total_minus, omitted_plus, omitted_minus
    real(dp) :: total_bar, omitted_bar, lhs, rhs, derivative_error, adjoint_error
    logical :: all_passed

    all_passed = .true.
    degree_indices = [0, 2, 4, 3]
    orders = [0, 1, 4, 2]
    coefficients = [ &
        cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.2_dp, 0.6_dp, dp), cmplx(-0.1_dp, 0.5_dp, dp)]
    coefficients_dot = [ &
        cmplx(0.03_dp, 0.01_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp), &
        cmplx(0.01_dp, -0.03_dp, dp), cmplx(0.02_dp, -0.01_dp, dp)]
    degree_limit = 3
    order_limit = 2

    call analyze_toroidal_spectral_modes( &
        degree_indices, orders, coefficients, degree_limit, order_limit, .true., &
        retained_count, omitted_count, zero_mode_count, total_energy, omitted_energy, &
        status)
    call record_condition(status == 0 .and. retained_count == 3 .and. &
        omitted_count == 1 .and. zero_mode_count == 1 .and. &
        abs(total_energy - sum(abs(coefficients)**2)) < 2.0e-14_dp .and. &
        abs(omitted_energy - abs(coefficients(3))**2) < 2.0e-14_dp, &
        "toroidal mode metadata reports rectangular truncation and zero mode")
    call analyze_toroidal_spectral_modes( &
        degree_indices, orders, coefficients, degree_limit, order_limit, .false., &
        retained_count, omitted_count, zero_mode_count, total_energy, omitted_energy, &
        status)
    call record_condition(status /= 0, &
        "toroidal mode metadata can reject an explicit zero mode")

    call analyze_toroidal_spectral_modes_jvp( &
        degree_indices, orders, coefficients, degree_limit, order_limit, .true., &
        coefficients_dot, total_energy_dot, omitted_energy_dot, status)
    call analyze_toroidal_spectral_modes( &
        degree_indices, orders, coefficients + step*coefficients_dot, degree_limit, &
        order_limit, .true., retained_count, omitted_count, zero_mode_count, total_plus, &
        omitted_plus, status)
    call analyze_toroidal_spectral_modes( &
        degree_indices, orders, coefficients - step*coefficients_dot, degree_limit, &
        order_limit, .true., retained_count, omitted_count, zero_mode_count, total_minus, &
        omitted_minus, status)
    derivative_error = max(abs(total_energy_dot - &
        (total_plus - total_minus)/(2.0_dp*step)), abs(omitted_energy_dot - &
        (omitted_plus - omitted_minus)/(2.0_dp*step)))
    call record_condition(status == 0 .and. derivative_error < 2.0e-8_dp, &
        "toroidal truncation energy JVP matches central reassembly")

    total_bar = 0.8_dp
    omitted_bar = -0.3_dp
    call analyze_toroidal_spectral_modes_vjp( &
        degree_indices, orders, coefficients, degree_limit, order_limit, .true., &
        total_bar, omitted_bar, coefficients_bar, status)
    lhs = total_bar*total_energy_dot + omitted_bar*omitted_energy_dot
    rhs = real(sum(conjg(coefficients_bar)*coefficients_dot), dp)
    adjoint_error = abs(lhs - rhs)
    call record_condition(status == 0 .and. adjoint_error < &
        5.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "toroidal truncation energy VJP satisfies the real complex adjoint")

    call check_summary("toroidal spectral metadata")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_toroidal_spectral_metadata_ad
