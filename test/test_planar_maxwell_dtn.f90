program test_planar_maxwell_dtn
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: apply_planar_maxwell_dtn
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: nx = 8, ny = 6
    real(dp), parameter :: length_x = 2.0_dp*acos(-1.0_dp)
    real(dp), parameter :: length_y = 3.0_dp
    real(dp), parameter :: wave_number = 4.2_dp
    complex(dp) :: derivative(2, nx, ny), expected(2, nx, ny)
    complex(dp) :: phase, trace(2, nx, ny)
    real(dp) :: beta, x
    integer :: i, j, status
    logical :: all_passed

    all_passed = .true.
    trace = cmplx(0.0_dp, 0.0_dp, dp)
    trace(1, :, :) = cmplx(1.0_dp, 0.0_dp, dp)
    call apply_planar_maxwell_dtn( &
        trace, wave_number, length_x, length_y, derivative, status)
    expected = cmplx(0.0_dp, 0.0_dp, dp)
    expected(1, :, :) = cmplx(0.0_dp, -wave_number, dp)
    call record_condition(status == 0 .and. &
        maxval(abs(derivative - expected)) < 2.0e-12_dp, &
        "normal-incidence Maxwell DtN has outgoing impedance")

    beta = sqrt(wave_number**2 - 1.0_dp)
    trace = cmplx(0.0_dp, 0.0_dp, dp)
    expected = cmplx(0.0_dp, 0.0_dp, dp)
    do j = 1, ny
        do i = 1, nx
            x = length_x*real(i - 1, dp)/real(nx, dp)
            phase = exp(cmplx(0.0_dp, x, dp))
            trace(2, i, j) = phase
            expected(2, i, j) = cmplx(0.0_dp, -beta, dp)*phase
        end do
    end do
    call apply_planar_maxwell_dtn( &
        trace, wave_number, length_x, length_y, derivative, status)
    call record_condition(maxval(abs(derivative - expected)) < 3.0e-12_dp, &
        "transverse electric Fourier mode has the TE eigenvalue")

    trace = cmplx(0.0_dp, 0.0_dp, dp)
    expected = cmplx(0.0_dp, 0.0_dp, dp)
    do j = 1, ny
        do i = 1, nx
            x = length_x*real(i - 1, dp)/real(nx, dp)
            phase = exp(cmplx(0.0_dp, x, dp))
            trace(1, i, j) = phase
            expected(1, i, j) = &
                cmplx(0.0_dp, -wave_number**2/beta, dp)*phase
        end do
    end do
    call apply_planar_maxwell_dtn( &
        trace, wave_number, length_x, length_y, derivative, status)
    call record_condition(maxval(abs(derivative - expected)) < 3.0e-12_dp, &
        "longitudinal electric Fourier mode has the TM eigenvalue")

    call check_summary("Planar Maxwell DtN")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_planar_maxwell_dtn
