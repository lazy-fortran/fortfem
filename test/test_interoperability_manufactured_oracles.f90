program test_interoperability_manufactured_oracles
    use check, only: check_condition, check_summary
    use fortfem_generated_interoperability_oracles, only: &
        generated_interoperability_oracles
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: values(7)
    logical :: all_passed

    all_passed = .true.
    call generated_interoperability_oracles(0.25_dp, 0.4_dp, values)
    call record_condition( &
        abs(values(1) - sin(0.25_dp*pi)*sin(0.4_dp*pi)) < 1.0e-13_dp, &
        "Generated Poisson exact field matches the analytical oracle")
    call record_condition( &
        abs(values(2) - 2.0_dp*pi**2*values(1)) < 1.0e-12_dp, &
        "Generated Poisson forcing is minus the exact Laplacian")
    call record_condition( &
        maxval(abs(values(3:4) - [sin(0.4_dp*pi), sin(0.25_dp*pi)])) < &
        1.0e-13_dp, &
        "Generated Ampere field matches the analytical oracle")
    call record_condition( &
        abs(values(5) - pi*(cos(0.25_dp*pi) - cos(0.4_dp*pi))) < &
        1.0e-12_dp, &
        "Generated Ampere scalar curl matches the analytical oracle")
    call record_condition( &
        maxval(abs(values(6:7) - (pi**2 + 1.0_dp)*values(3:4))) < &
        1.0e-11_dp, &
        "Generated Ampere forcing equals curl curl E plus E")

    call check_zero_tangential_trace(0.37_dp)
    call check_summary("Interoperability manufactured analytical oracles")
    if (.not. all_passed) error stop 1

contains

    subroutine check_zero_tangential_trace(parameter)
        real(dp), intent(in) :: parameter

        real(dp) :: boundary_values(7)

        call generated_interoperability_oracles(0.0_dp, parameter, boundary_values)
        call record_condition(abs(boundary_values(4)) < 1.0e-14_dp, &
            "Ampere tangential trace vanishes on x=0")
        call generated_interoperability_oracles(1.0_dp, parameter, boundary_values)
        call record_condition(abs(boundary_values(4)) < 1.0e-14_dp, &
            "Ampere tangential trace vanishes on x=1")
        call generated_interoperability_oracles(parameter, 0.0_dp, boundary_values)
        call record_condition(abs(boundary_values(3)) < 1.0e-14_dp, &
            "Ampere tangential trace vanishes on y=0")
        call generated_interoperability_oracles(parameter, 1.0_dp, boundary_values)
        call record_condition(abs(boundary_values(3)) < 1.0e-14_dp, &
            "Ampere tangential trace vanishes on y=1")
    end subroutine check_zero_tangential_trace

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_interoperability_manufactured_oracles
