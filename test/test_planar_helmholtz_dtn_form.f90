program test_planar_helmholtz_dtn_form
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_planar_helmholtz_dtn_form
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: sample_count = 11
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: period = 2.0_dp*pi
    complex(dp) :: form(sample_count, sample_count)
    complex(dp) :: mode_values(sample_count), result(sample_count)
    complex(dp) :: expected_eigenvalue
    real(dp) :: mass_eigenvalue
    integer :: mode, point, status
    logical :: all_passed

    all_passed = .true.
    call assemble_planar_helmholtz_dtn_form( &
        sample_count, 3.0_dp, period, form, status)
    call record_condition(status == 0, "Planar DtN form assembly succeeds")

    do mode = 0, 4
        do point = 1, sample_count
            mode_values(point) = exp(cmplx(0.0_dp, &
                2.0_dp*pi*real(mode*(point - 1), dp)/ &
                real(sample_count, dp), dp))
        end do
        if (mode <= 3) then
            expected_eigenvalue = cmplx( &
                0.0_dp, sqrt(9.0_dp - real(mode**2, dp)), dp)
        else
            expected_eigenvalue = cmplx( &
                -sqrt(real(mode**2 - 9, dp)), 0.0_dp, dp)
        end if
        mass_eigenvalue = period/real(sample_count, dp)/3.0_dp*( &
            2.0_dp + cos(2.0_dp*pi*real(mode, dp)/ &
            real(sample_count, dp)))
        result = matmul(form, mode_values)
        call record_condition(maxval(abs(result - &
            mass_eigenvalue*expected_eigenvalue*mode_values)) < 2.0e-12_dp, &
            "DtN weak form has the analytical Fourier eigenvalue")
    end do

    call check_summary("Planar Helmholtz DtN weak form")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition
end program test_planar_helmholtz_dtn_form
