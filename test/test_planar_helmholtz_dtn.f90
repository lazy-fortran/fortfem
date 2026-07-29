program test_planar_helmholtz_dtn
    use fortfem_kinds, only: dp
    use fortfem_planar_helmholtz_dtn, only: apply_planar_helmholtz_dtn
    use check, only: check_condition, check_summary
    implicit none

    logical :: all_passed

    all_passed = .true.
    call test_fourier_mode(0, 2.0_dp, cmplx(0.0_dp, 2.0_dp, dp), &
        "Constant propagating mode has outgoing imaginary derivative")
    call test_fourier_mode(1, 2.0_dp, &
        cmplx(0.0_dp, sqrt(3.0_dp), dp), &
        "Nonzero propagating mode has outgoing imaginary derivative")
    call test_fourier_mode(3, 2.0_dp, &
        cmplx(-sqrt(5.0_dp), 0.0_dp, dp), &
        "Evanescent mode decays into the exterior")
    call test_fourier_mode(-3, 2.0_dp, &
        cmplx(-sqrt(5.0_dp), 0.0_dp, dp), &
        "Negative evanescent mode uses the same outgoing branch")
    call check_summary("Planar Helmholtz DtN")
    if (.not. all_passed) error stop 1

contains

    subroutine test_fourier_mode(mode, wavenumber, eigenvalue, description)
        integer, intent(in) :: mode
        real(dp), intent(in) :: wavenumber
        complex(dp), intent(in) :: eigenvalue
        character(len=*), intent(in) :: description

        integer, parameter :: n = 7
        real(dp), parameter :: period = 2.0_dp * acos(-1.0_dp)
        complex(dp) :: trace(n), derivative(n)
        real(dp) :: coordinate
        integer :: point_id

        do point_id = 1, n
            coordinate = period * real(point_id - 1, dp) / real(n, dp)
            trace(point_id) = exp(cmplx(0.0_dp, &
                real(mode, dp) * coordinate, dp))
        end do

        call apply_planar_helmholtz_dtn( &
            trace, wavenumber, period, derivative)

        call record_condition( &
            maxval(abs(derivative - eigenvalue * trace)) < 1.0e-12_dp, &
            description)
    end subroutine test_fourier_mode

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_planar_helmholtz_dtn
