program test_toroidal_poisson_dtn
    use check, only: check_condition, check_summary
    use fortfem_fourier, only: &
        evaluate_toroidal_harmonic_p, toroidal_poisson_exterior_dtn_p
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: scale = 1.7_dp
    real(dp), parameter :: eta = &
        1.3169578969248167086250463473079684_dp
    real(dp), parameter :: theta = 0.6_dp
    real(dp), parameter :: phi = 0.4_dp
    real(dp), parameter :: reference_value = &
        1.6160590342216749833764963336096247_dp
    real(dp), parameter :: reference_normal = &
        -2.5743996635144626758125703168794369_dp
    real(dp), parameter :: reference_dtn = &
        -1.5930109043042124475276790110593581_dp
    real(dp) :: value, normal_derivative, dtn_value
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call evaluate_toroidal_harmonic_p( &
        2, 1, eta, theta, phi, value, status)
    call record(status == 0 .and. &
        abs(value - reference_value) < 5.0e-12_dp, &
        "Toroidal harmonic matches a 50-digit mpmath oracle")

    call toroidal_poisson_exterior_dtn_p( &
        2, 1, scale, eta, theta, phi, &
        value, normal_derivative, dtn_value, status)
    call record(status == 0 .and. &
        abs(normal_derivative - reference_normal) < 8.0e-12_dp, &
        "Toroidal Poisson normal derivative matches the oracle")
    call record(abs(dtn_value - reference_dtn) < 8.0e-12_dp, &
        "Toroidal Poisson DtN ratio matches the oracle")

    call toroidal_poisson_exterior_dtn_p( &
        2, 1, 0.0_dp, eta, theta, phi, &
        value, normal_derivative, dtn_value, status)
    call record(status /= 0, "Toroidal DtN rejects a zero scale")

    call check_summary("Toroidal Poisson DtN")
    if (.not. all_passed) error stop 1

contains

    subroutine record(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        call check_condition(condition, message)
        all_passed = all_passed .and. condition
    end subroutine record

end program test_toroidal_poisson_dtn
