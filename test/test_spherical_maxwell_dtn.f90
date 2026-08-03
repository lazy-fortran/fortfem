program test_spherical_maxwell_dtn
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        apply_spherical_maxwell_dtn, spherical_maxwell_dtn_eigenvalues
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: k = 2.0_dp
    real(dp), parameter :: radius = 1.5_dp
    complex(dp), parameter :: imaginary_unit = cmplx(0.0_dp, 1.0_dp, dp)
    complex(dp) :: logarithmic_derivative, te_eigenvalue, tm_eigenvalue
    complex(dp) :: te_trace(2), tm_trace(2), te_result(2), tm_result(2)
    integer :: degrees(2), status
    logical :: all_passed

    all_passed = .true.
    call spherical_maxwell_dtn_eigenvalues( &
        1, k, radius, te_eigenvalue, tm_eigenvalue, status)
    logarithmic_derivative = k*( &
        imaginary_unit + 1.0_dp/(k*radius + imaginary_unit) - &
        2.0_dp/(k*radius))
    call record_condition(status == 0 .and. abs(te_eigenvalue + &
        logarithmic_derivative + 1.0_dp/radius) < 2.0e-14_dp, &
        "Spherical Maxwell TE mode matches the Riccati-Hankel derivative")
    call record_condition(abs( &
        te_eigenvalue*tm_eigenvalue + k**2) < 2.0e-14_dp, &
        "Spherical Maxwell TE and TM impedances satisfy exact duality")
    call record_condition(aimag(te_eigenvalue) < 0.0_dp .and. &
        aimag(tm_eigenvalue) < 0.0_dp, &
        "Spherical Maxwell DtN selects outgoing passive modes")

    degrees = [1, 2]
    te_trace = [cmplx(1.0_dp, -0.5_dp, dp), cmplx(0.25_dp, 0.75_dp, dp)]
    tm_trace = [cmplx(-0.2_dp, 0.4_dp, dp), cmplx(0.5_dp, -0.1_dp, dp)]
    call apply_spherical_maxwell_dtn( &
        degrees, te_trace, tm_trace, k, radius, te_result, tm_result, status)
    call record_condition(status == 0 .and. &
        abs(te_result(1)/te_trace(1) - te_eigenvalue) < 2.0e-14_dp .and. &
        abs(tm_result(1)/tm_trace(1) - tm_eigenvalue) < 2.0e-14_dp, &
        "Spherical Maxwell DtN applies diagonal TE and TM modal maps")

    call spherical_maxwell_dtn_eigenvalues( &
        0, k, radius, te_eigenvalue, tm_eigenvalue, status)
    call record_condition(status /= 0, &
        "Spherical Maxwell DtN rejects the nonexistent degree-zero mode")
    call check_summary("Spherical Maxwell DtN")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_spherical_maxwell_dtn
