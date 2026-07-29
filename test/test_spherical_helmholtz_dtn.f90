program test_spherical_helmholtz_dtn
    use fortfem_api, only: apply_spherical_helmholtz_dtn, &
        spherical_helmholtz_dtn_eigenvalue
    use fortfem_kinds, only: dp
    use check, only: check_condition, check_summary
    implicit none

    real(dp), parameter :: k = 2.0_dp
    real(dp), parameter :: radius = 1.5_dp
    complex(dp), parameter :: imaginary_unit = cmplx(0.0_dp, 1.0_dp, dp)
    complex(dp) :: eigenvalue, expected
    complex(dp) :: trace(3), normal_derivative(3)
    integer :: degrees(3), status
    logical :: all_passed

    all_passed = .true.
    call spherical_helmholtz_dtn_eigenvalue(0, k, radius, eigenvalue, status)
    expected = imaginary_unit*k - 1.0_dp/radius
    call record_condition(status == 0 .and. &
        abs(eigenvalue - expected) < 2.0e-14_dp, &
        "Spherical degree-zero DtN matches exp(ikr)/r")

    call spherical_helmholtz_dtn_eigenvalue(1, k, radius, eigenvalue, status)
    expected = k*(imaginary_unit + 1.0_dp/(k*radius + imaginary_unit) - &
        2.0_dp/(k*radius))
    call record_condition(status == 0 .and. &
        abs(eigenvalue - expected) < 2.0e-14_dp, &
        "Spherical degree-one DtN matches the closed Hankel form")

    degrees = [0, 1, 0]
    trace = [cmplx(1.0_dp, -0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp), &
        cmplx(0.5_dp, 0.125_dp, dp)]
    call apply_spherical_helmholtz_dtn( &
        degrees, trace, k, radius, normal_derivative, status)
    call record_condition( &
        status == 0, "Spherical DtN applies to modal coefficients")
    call record_condition(abs(normal_derivative(1)/trace(1) - &
        normal_derivative(3)/trace(3)) < 2.0e-14_dp, &
        "Spherical DtN is degenerate across equal degrees")

    call spherical_helmholtz_dtn_eigenvalue(-1, k, radius, &
        eigenvalue, status)
    call record_condition(status /= 0, "Spherical DtN rejects negative degree")

    call check_summary("Spherical Helmholtz DtN")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        call check_condition(condition, description)
        all_passed = all_passed .and. condition
    end subroutine record_condition
end program test_spherical_helmholtz_dtn
