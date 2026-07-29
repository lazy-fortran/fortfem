program test_cartesian_helmholtz_pml
    use check, only: check_condition, check_summary
    use fortfem_api, only: cartesian_curl_curl_pml_coefficients, &
        cartesian_scalar_helmholtz_pml_coefficients
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), parameter :: one = cmplx(1.0_dp, 0.0_dp, dp)
    complex(dp), parameter :: stretch = cmplx(1.0_dp, 2.0_dp, dp)
    complex(dp) :: curl_coefficient(3), gradient_coefficient(3)
    complex(dp) :: mass, mass_coefficient(3), stretches(3)
    logical :: all_passed
    integer :: status

    all_passed = .true.
    stretches = [stretch, one, one]

    call cartesian_scalar_helmholtz_pml_coefficients( &
        stretches, gradient_coefficient, mass, status)
    call record_condition(status == 0, "Scalar Helmholtz PML accepts x stretch")
    call record_condition(maxval(abs(gradient_coefficient - &
        [one/stretch, stretch, stretch])) < 2.0e-15_dp, &
        "Scalar PML gradient tensor matches the coordinate transform")
    call record_condition(abs(mass - stretch) < 2.0e-15_dp, &
        "Scalar PML mass Jacobian matches the coordinate transform")

    call cartesian_curl_curl_pml_coefficients( &
        stretches, curl_coefficient, mass_coefficient, status)
    call record_condition(status == 0, "Curl-curl PML accepts x stretch")
    call record_condition(maxval(abs(curl_coefficient - &
        [stretch, one/stretch, one/stretch])) < 2.0e-15_dp, &
        "Curl-curl PML tensor follows the covariant Piola transform")
    call record_condition(maxval(abs(mass_coefficient - &
        [one/stretch, stretch, stretch])) < 2.0e-15_dp, &
        "Maxwell PML mass tensor follows the covariant Piola transform")

    stretches(2) = cmplx(0.0_dp, 0.0_dp, dp)
    call cartesian_scalar_helmholtz_pml_coefficients( &
        stretches, gradient_coefficient, mass, status)
    call record_condition(status /= 0, "PML rejects a singular stretch")

    call check_summary("Cartesian Helmholtz PML coefficients")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        call check_condition(condition, description)
        all_passed = all_passed .and. condition
    end subroutine record_condition
end program test_cartesian_helmholtz_pml
