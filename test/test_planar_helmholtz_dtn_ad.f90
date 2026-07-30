program test_planar_helmholtz_dtn_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: apply_planar_helmholtz_dtn, &
        apply_planar_helmholtz_dtn_jvp, apply_planar_helmholtz_dtn_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: sample_count = 9
    real(dp), parameter :: difference_step = 1.0e-6_dp
    real(dp), parameter :: period = 2.0_dp*acos(-1.0_dp)
    real(dp), parameter :: wavenumber = 3.3_dp
    complex(dp) :: derivative_bar(sample_count)
    complex(dp) :: derivative_dot(sample_count)
    complex(dp) :: derivative_minus(sample_count), derivative_plus(sample_count)
    complex(dp) :: trace(sample_count), trace_bar(sample_count)
    complex(dp) :: trace_dot(sample_count)
    real(dp) :: forward_pairing, period_bar, period_dot
    real(dp) :: reverse_pairing, wavenumber_bar, wavenumber_dot
    integer :: point
    logical :: all_passed

    all_passed = .true.
    do point = 1, sample_count
        trace(point) = cmplx( &
            sin(0.37_dp*real(point, dp)), cos(0.23_dp*real(point, dp)), dp)
        trace_dot(point) = cmplx( &
            0.1_dp*cos(0.41_dp*real(point, dp)), &
            -0.07_dp*sin(0.29_dp*real(point, dp)), dp)
        derivative_bar(point) = cmplx( &
            -0.08_dp*sin(0.31_dp*real(point, dp)), &
            0.09_dp*cos(0.19_dp*real(point, dp)), dp)
    end do
    wavenumber_dot = 0.06_dp
    period_dot = -0.04_dp

    call apply_planar_helmholtz_dtn_jvp( &
        trace, wavenumber, period, trace_dot, wavenumber_dot, period_dot, &
        derivative_dot)
    call apply_planar_helmholtz_dtn( &
        trace - difference_step*trace_dot, &
        wavenumber - difference_step*wavenumber_dot, &
        period - difference_step*period_dot, derivative_minus)
    call apply_planar_helmholtz_dtn( &
        trace + difference_step*trace_dot, &
        wavenumber + difference_step*wavenumber_dot, &
        period + difference_step*period_dot, derivative_plus)
    call record_condition(maxval(abs(derivative_dot - &
        (derivative_plus - derivative_minus)/(2.0_dp*difference_step))) < &
        3.0e-9_dp, "Planar Helmholtz DtN JVP matches central difference")

    call apply_planar_helmholtz_dtn_vjp( &
        trace, wavenumber, period, derivative_bar, trace_bar, &
        wavenumber_bar, period_bar)
    forward_pairing = real(sum(conjg(derivative_bar)*derivative_dot), dp)
    reverse_pairing = real(sum(conjg(trace_bar)*trace_dot), dp) + &
        wavenumber_bar*wavenumber_dot + period_bar*period_dot
    call record_condition(abs(forward_pairing - reverse_pairing) < 3.0e-12_dp, &
        "Planar Helmholtz DtN JVP and VJP satisfy the complex dot identity")

    call check_summary("Planar Helmholtz DtN derivatives")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_planar_helmholtz_dtn_ad
