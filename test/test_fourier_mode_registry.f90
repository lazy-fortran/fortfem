program test_fourier_mode_registry
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_fourier_mode, evaluate_fourier_mode_jvp, &
        evaluate_fourier_mode_vjp, find_fourier_mode, &
        fourier_mode_conjugate_index, fourier_mode_triad_closed, &
        fourier_mode_registry_t, initialize_fourier_mode_registry, &
        validate_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: mode_count = 5
    integer, parameter :: poloidal(mode_count) = [-1, 0, 1, 2, -2]
    integer, parameter :: toroidal(mode_count) = [1, 0, -1, 2, -2]
    integer, parameter :: radial_power(mode_count) = [1, 0, 1, 2, 2]
    real(dp), parameter :: normalization(mode_count) = [ &
        1.0_dp, 2.0_dp, 1.0_dp, 0.5_dp, 0.5_dp]
    integer, parameter :: field_periods = 3
    real(dp), parameter :: poloidal_phase = 0.1_dp
    real(dp), parameter :: toroidal_phase = -0.2_dp
    real(dp), parameter :: radius = 0.7_dp
    real(dp), parameter :: theta = 0.4_dp
    real(dp), parameter :: phi = -0.3_dp
    real(dp), parameter :: radius_dot = 0.13_dp
    real(dp), parameter :: theta_dot = -0.17_dp
    real(dp), parameter :: phi_dot = 0.11_dp
    real(dp), parameter :: eps = 1.0e-7_dp
    complex(dp), parameter :: value_bar = cmplx(0.7_dp, -0.4_dp, dp)

    type(fourier_mode_registry_t) :: registry, copied
    type(fortsparse_status_t) :: status
    complex(dp) :: value, value_dot, value_plus, value_minus
    complex(dp) :: expected, expected_dot
    complex(dp) :: radial_derivative, theta_derivative, phi_derivative
    real(dp) :: radius_bar, theta_bar, phi_bar, lhs, rhs
    real(dp) :: expected_radius_bar, expected_theta_bar, expected_phi_bar
    integer :: mode, conjugate, triad_output

    call initialize_fourier_mode_registry( &
        registry, poloidal, toroidal, field_periods, poloidal_phase, &
        toroidal_phase, .true., radial_power, normalization, status)
    call check_condition(status%code == 0, &
        "Fourier mode registry accepts paired real-packed modes")
    call check_condition(validate_fourier_mode_registry(registry, status) .and. &
        status%code == 0, "Fourier mode registry validates its metadata")

    copied = registry
    copied%poloidal_modes(1) = 99
    call check_condition(registry%poloidal_modes(1) == poloidal(1) .and. &
        copied%poloidal_modes(1) == 99, &
        "Fourier mode registry assignment is a deep copy")

    mode = find_fourier_mode(registry, -1, 1)
    conjugate = fourier_mode_conjugate_index(registry, mode, status)
    call check_condition(status%code == 0 .and. mode == 1 .and. conjugate == 3, &
        "Fourier registry finds conjugate mode pairs")
    triad_output = find_fourier_mode(registry, 0, 0)
    call check_condition(fourier_mode_triad_closed( &
        registry, mode, conjugate, triad_output, status) .and. status%code == 0, &
        "Fourier registry accepts a retained triad")
    call check_condition(.not. fourier_mode_triad_closed( &
        registry, mode, mode, triad_output, status) .and. status%code == 0, &
        "Fourier registry reports an absent triad")

    mode = find_fourier_mode(registry, 1, -1)
    call evaluate_fourier_mode( &
        registry, mode, radius, theta, phi, value, radial_derivative, &
        theta_derivative, phi_derivative, status)
    expected = normalization(mode)*radius**radial_power(mode)*exp(cmplx( &
        0.0_dp, real(poloidal(mode), dp)*(theta + poloidal_phase) + &
        real(toroidal(mode), dp)*real(field_periods, dp)* &
        (phi + toroidal_phase), dp))
    call check_condition(status%code == 0 .and. &
        abs(value - expected) < 1.0e-14_dp, &
        "Fourier registry evaluates phase and radial factors")

    call evaluate_fourier_mode_jvp( &
        registry, mode, radius, theta, phi, radius_dot, theta_dot, phi_dot, &
        value_dot, status)
    expected_dot = radial_derivative*radius_dot + theta_derivative*theta_dot + &
        phi_derivative*phi_dot
    call evaluate_fourier_mode( &
        registry, mode, radius + eps*radius_dot, theta + eps*theta_dot, &
        phi + eps*phi_dot, value_plus, radial_derivative, theta_derivative, &
        phi_derivative, status)
    call evaluate_fourier_mode( &
        registry, mode, radius - eps*radius_dot, theta - eps*theta_dot, &
        phi - eps*phi_dot, value_minus, radial_derivative, theta_derivative, &
        phi_derivative, status)
    call check_condition(status%code == 0 .and. &
        abs(value_dot - expected_dot) < 1.0e-14_dp .and. &
        abs((value_plus - value_minus)/(2.0_dp*eps) - value_dot) < 1.0e-8_dp, &
        "Fourier registry JVP matches central differences")

    call evaluate_fourier_mode_vjp( &
        registry, mode, radius, theta, phi, value_bar, radius_bar, theta_bar, &
        phi_bar, status)
    call evaluate_fourier_mode( &
        registry, mode, radius, theta, phi, value, radial_derivative, &
        theta_derivative, phi_derivative, status)
    lhs = real(conjg(value_bar)*value_dot, dp)
    rhs = radius_bar*radius_dot + theta_bar*theta_dot + phi_bar*phi_dot
    expected_radius_bar = real(conjg(value_bar)*radial_derivative, dp)
    expected_theta_bar = real(conjg(value_bar)*theta_derivative, dp)
    expected_phi_bar = real(conjg(value_bar)*phi_derivative, dp)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp .and. &
        abs(radius_bar - expected_radius_bar) < 1.0e-13_dp .and. &
        abs(theta_bar - expected_theta_bar) < 1.0e-13_dp .and. &
        abs(phi_bar - expected_phi_bar) < 1.0e-13_dp, &
        "Fourier registry VJP satisfies the complex real-part identity")

    call initialize_fourier_mode_registry( &
        registry, [0, 1], [0, 2], field_periods, 0.0_dp, 0.0_dp, .true., &
        status=status)
    call check_condition(status%code /= 0, &
        "Fourier registry rejects unpaired real-packed modes")
    call initialize_fourier_mode_registry( &
        registry, [0, 0], [0, 0], field_periods, 0.0_dp, 0.0_dp, .false., &
        status=status)
    call check_condition(status%code /= 0, &
        "Fourier registry rejects duplicate modes")
    call initialize_fourier_mode_registry( &
        registry, [0, 1], [0, 1], field_periods, 0.0_dp, 0.0_dp, .false., &
        radial_powers=[0], status=status)
    call check_condition(status%code /= 0, &
        "Fourier registry rejects a mismatched radial-power array")

    call evaluate_fourier_mode( &
        copied, 0, radius, theta, phi, value, radial_derivative, &
        theta_derivative, phi_derivative, status)
    call check_condition(status%code /= 0 .and. &
        ieee_is_finite(real(value, dp)), &
        "Fourier registry rejects an out-of-range mode index")

    call check_summary("Fourier mode registry")
end program test_fourier_mode_registry
