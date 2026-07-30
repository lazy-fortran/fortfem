program test_cartesian_helmholtz_pml
    use check, only: check_condition, check_summary
    use fortfem_api, only: cartesian_curl_curl_pml_coefficients, &
        cartesian_curl_curl_pml_coefficients_jvp, &
        cartesian_curl_curl_pml_coefficients_vjp, &
        cartesian_scalar_helmholtz_pml_coefficients, &
        cartesian_scalar_helmholtz_pml_coefficients_jvp, &
        cartesian_scalar_helmholtz_pml_coefficients_vjp
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), parameter :: one = cmplx(1.0_dp, 0.0_dp, dp)
    complex(dp), parameter :: stretch = cmplx(1.0_dp, 2.0_dp, dp)
    complex(dp) :: curl_coefficient(3), gradient_coefficient(3)
    complex(dp) :: curl_bar(3), curl_dot(3), curl_minus(3), curl_plus(3)
    complex(dp) :: gradient_bar(3), gradient_dot(3)
    complex(dp) :: gradient_minus(3), gradient_plus(3)
    complex(dp) :: mass, mass_coefficient(3), stretches(3)
    complex(dp) :: mass_bar, mass_coefficient_bar(3)
    complex(dp) :: mass_dot, mass_minus, mass_plus
    complex(dp) :: mass_coefficient_dot(3)
    complex(dp) :: mass_coefficient_minus(3), mass_coefficient_plus(3)
    complex(dp) :: stretch_bar(3), stretch_dot(3)
    logical :: all_passed
    real(dp), parameter :: difference_step = 1.0e-6_dp
    real(dp) :: forward_pairing, reverse_pairing
    integer :: status, status_minus, status_plus

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

    stretches = [ &
        cmplx(1.1_dp, 0.3_dp, dp), cmplx(0.9_dp, 0.2_dp, dp), &
        cmplx(1.2_dp, 0.1_dp, dp)]
    stretch_dot = [ &
        cmplx(0.07_dp, -0.02_dp, dp), cmplx(-0.03_dp, 0.04_dp, dp), &
        cmplx(0.02_dp, 0.01_dp, dp)]
    gradient_bar = [ &
        cmplx(0.4_dp, -0.1_dp, dp), cmplx(-0.2_dp, 0.3_dp, dp), &
        cmplx(0.1_dp, 0.5_dp, dp)]
    mass_bar = cmplx(-0.3_dp, 0.2_dp, dp)
    call cartesian_scalar_helmholtz_pml_coefficients_jvp( &
        stretches, stretch_dot, gradient_dot, mass_dot, status)
    call cartesian_scalar_helmholtz_pml_coefficients( &
        stretches - difference_step*stretch_dot, gradient_minus, mass_minus, &
        status_minus)
    call cartesian_scalar_helmholtz_pml_coefficients( &
        stretches + difference_step*stretch_dot, gradient_plus, mass_plus, &
        status_plus)
    call record_condition(status == 0 .and. status_minus == 0 .and. &
        status_plus == 0, "Scalar PML JVP accepts nonsingular stretches")
    call record_condition(maxval(abs(gradient_dot - &
        (gradient_plus - gradient_minus)/(2.0_dp*difference_step))) < &
        2.0e-9_dp, "Scalar PML gradient JVP matches central difference")
    call record_condition(abs(mass_dot - &
        (mass_plus - mass_minus)/(2.0_dp*difference_step)) < 2.0e-9_dp, &
        "Scalar PML mass JVP matches central difference")
    call cartesian_scalar_helmholtz_pml_coefficients_vjp( &
        stretches, gradient_bar, mass_bar, stretch_bar, status)
    forward_pairing = real(sum(conjg(gradient_bar)*gradient_dot) + &
        conjg(mass_bar)*mass_dot, dp)
    reverse_pairing = real(sum(conjg(stretch_bar)*stretch_dot), dp)
    call record_condition(abs(forward_pairing - reverse_pairing) < 2.0e-12_dp, &
        "Scalar PML JVP and VJP satisfy the real complex dot identity")

    curl_bar = [ &
        cmplx(-0.1_dp, 0.2_dp, dp), cmplx(0.3_dp, 0.1_dp, dp), &
        cmplx(0.2_dp, -0.4_dp, dp)]
    mass_coefficient_bar = [ &
        cmplx(0.5_dp, 0.2_dp, dp), cmplx(-0.1_dp, 0.3_dp, dp), &
        cmplx(0.2_dp, 0.1_dp, dp)]
    call cartesian_curl_curl_pml_coefficients_jvp( &
        stretches, stretch_dot, curl_dot, mass_coefficient_dot, status)
    call cartesian_curl_curl_pml_coefficients( &
        stretches - difference_step*stretch_dot, curl_minus, &
        mass_coefficient_minus, status_minus)
    call cartesian_curl_curl_pml_coefficients( &
        stretches + difference_step*stretch_dot, curl_plus, &
        mass_coefficient_plus, status_plus)
    call record_condition(maxval(abs(curl_dot - &
        (curl_plus - curl_minus)/(2.0_dp*difference_step))) < 2.0e-9_dp, &
        "Curl-curl PML JVP matches central difference")
    call record_condition(maxval(abs(mass_coefficient_dot - &
        (mass_coefficient_plus - mass_coefficient_minus)/ &
        (2.0_dp*difference_step))) < 2.0e-9_dp, &
        "Maxwell PML mass JVP matches central difference")
    call cartesian_curl_curl_pml_coefficients_vjp( &
        stretches, curl_bar, mass_coefficient_bar, stretch_bar, status)
    forward_pairing = real(sum(conjg(curl_bar)*curl_dot) + &
        sum(conjg(mass_coefficient_bar)*mass_coefficient_dot), dp)
    reverse_pairing = real(sum(conjg(stretch_bar)*stretch_dot), dp)
    call record_condition(abs(forward_pairing - reverse_pairing) < 2.0e-12_dp, &
        "Curl-curl PML JVP and VJP satisfy the real complex dot identity")

    stretches(2) = cmplx(0.0_dp, 0.0_dp, dp)
    call cartesian_scalar_helmholtz_pml_coefficients( &
        stretches, gradient_coefficient, mass, status)
    call record_condition(status /= 0, "PML rejects a singular stretch")
    call cartesian_scalar_helmholtz_pml_coefficients_jvp( &
        stretches, stretch_dot, gradient_dot, mass_dot, status)
    call record_condition(status /= 0 .and. &
        maxval(abs(gradient_dot)) == 0.0_dp .and. abs(mass_dot) == 0.0_dp, &
        "Scalar PML JVP rejects a singular stretch with zero outputs")
    call cartesian_curl_curl_pml_coefficients_vjp( &
        stretches, curl_bar, mass_coefficient_bar, stretch_bar, status)
    call record_condition(status /= 0 .and. maxval(abs(stretch_bar)) == 0.0_dp, &
        "Curl-curl PML VJP rejects a singular stretch with zero outputs")

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
