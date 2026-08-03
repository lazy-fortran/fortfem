program test_planar_helmholtz_dtn_form_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: assemble_planar_helmholtz_dtn_form, &
        assemble_planar_helmholtz_dtn_form_jvp, &
        assemble_planar_helmholtz_dtn_form_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: sample_count = 9
    real(dp), parameter :: difference_step = 1.0e-6_dp
    real(dp), parameter :: period = 2.0_dp*acos(-1.0_dp)
    real(dp), parameter :: wavenumber = 3.3_dp
    complex(dp) :: form_bar(sample_count, sample_count)
    complex(dp) :: form_dot(sample_count, sample_count)
    complex(dp) :: form_minus(sample_count, sample_count)
    complex(dp) :: form_plus(sample_count, sample_count)
    real(dp) :: forward_pairing, period_bar, period_dot
    real(dp) :: reverse_pairing, wavenumber_bar, wavenumber_dot
    integer :: column, row, status, status_minus, status_plus
    logical :: all_passed

    all_passed = .true.
    do column = 1, sample_count
        do row = 1, sample_count
            form_bar(row, column) = cmplx( &
                0.03_dp*sin(0.2_dp*real(row + 2*column, dp)), &
                -0.04_dp*cos(0.17_dp*real(2*row - column, dp)), dp)
        end do
    end do
    wavenumber_dot = 0.06_dp
    period_dot = -0.04_dp

    call assemble_planar_helmholtz_dtn_form_jvp( &
        sample_count, wavenumber, period, wavenumber_dot, period_dot, &
        form_dot, status)
    call assemble_planar_helmholtz_dtn_form( &
        sample_count, wavenumber - difference_step*wavenumber_dot, &
        period - difference_step*period_dot, form_minus, status_minus)
    call assemble_planar_helmholtz_dtn_form( &
        sample_count, wavenumber + difference_step*wavenumber_dot, &
        period + difference_step*period_dot, form_plus, status_plus)
    call record_condition(status == 0 .and. status_minus == 0 .and. &
        status_plus == 0, "Planar DtN form JVP accepts valid inputs")
    call record_condition(maxval(abs(form_dot - &
        (form_plus - form_minus)/(2.0_dp*difference_step))) < 3.0e-9_dp, &
        "Planar DtN form JVP matches central difference")

    call assemble_planar_helmholtz_dtn_form_vjp( &
        sample_count, wavenumber, period, form_bar, wavenumber_bar, &
        period_bar, status)
    forward_pairing = real(sum(conjg(form_bar)*form_dot), dp)
    reverse_pairing = &
        wavenumber_bar*wavenumber_dot + period_bar*period_dot
    call record_condition(status == 0, "Planar DtN form VJP succeeds")
    call record_condition(abs(forward_pairing - reverse_pairing) < 3.0e-12_dp, &
        "Planar DtN form JVP and VJP satisfy the complex dot identity")

    call check_summary("Planar Helmholtz DtN form derivatives")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_planar_helmholtz_dtn_form_ad
