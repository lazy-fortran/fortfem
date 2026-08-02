program test_fourier_facade
    !! Independent oracle for the canonical Fourier facade.
    !!
    !! The single retained mode is evaluated against its defining expression
    !! rather than against another FortFEM wrapper.  The toroidal P value is
    !! also checked against the pinned FortNum analytical reference.
    use check, only: check_condition, check_summary
    use fortfem_fourier, only: &
        evaluate_fourier_mode_expansion, fourier_mode_registry_t, &
        initialize_fourier_mode_registry, toroidal_p
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    integer, parameter :: mode_count = 1
    integer, parameter :: poloidal(mode_count) = [2]
    integer, parameter :: toroidal(mode_count) = [-1]
    integer, parameter :: radial_power(mode_count) = [3]
    real(dp), parameter :: normalization(mode_count) = [1.25_dp]
    integer, parameter :: field_periods = 3
    real(dp), parameter :: radius = 0.7_dp
    real(dp), parameter :: theta = -0.23_dp
    real(dp), parameter :: phi = 0.41_dp
    real(dp), parameter :: poloidal_phase = 0.17_dp
    real(dp), parameter :: toroidal_phase = -0.08_dp
    complex(dp), parameter :: coefficient = cmplx(0.6_dp, -0.9_dp, dp)
    complex(dp) :: coefficients(mode_count)
    complex(dp) :: value, radial_derivative, theta_derivative, phi_derivative
    complex(dp) :: expected, phase_factor
    type(fourier_mode_registry_t) :: registry
    type(fortsparse_status_t) :: status
    real(dp) :: phase, radial_factor

    coefficients = [coefficient]
    call initialize_fourier_mode_registry( &
        registry, poloidal, toroidal, field_periods, poloidal_phase, &
        toroidal_phase, .false., radial_powers=radial_power, &
        normalization=normalization, status=status)
    call check_condition(status%code == FORTSPARSE_OK, &
        "Fourier facade initializes a mode registry")

    call evaluate_fourier_mode_expansion( &
        registry, coefficients, radius, theta, phi, value, radial_derivative, &
        theta_derivative, phi_derivative, status)
    call check_condition(status%code == FORTSPARSE_OK, &
        "Fourier facade evaluates a single retained mode")

    phase = real(poloidal(1), dp)*(theta + poloidal_phase) + &
        real(toroidal(1), dp)*field_periods*(phi + toroidal_phase)
    phase_factor = exp(cmplx(0.0_dp, phase, dp))
    radial_factor = normalization(1)*radius**radial_power(1)
    expected = coefficient*radial_factor*phase_factor
    call check_condition(abs(value - expected) < 2.0e-13_dp, &
        "Fourier facade value agrees with the analytical mode expression")
    call check_condition(abs(radial_derivative - &
        coefficient*normalization(1)*real(radial_power(1), dp)* &
        radius**(radial_power(1) - 1)*phase_factor) < 2.0e-13_dp, &
        "Fourier facade radial derivative agrees with the power rule")
    call check_condition(abs(theta_derivative - &
        cmplx(0.0_dp, real(poloidal(1), dp), dp)*expected) < 2.0e-13_dp .and. &
        abs(phi_derivative - cmplx(0.0_dp, &
        real(toroidal(1), dp)*field_periods, dp)*expected) < 2.0e-13_dp, &
        "Fourier facade angular derivatives agree with the phase rule")

    call check_condition(abs(toroidal_p(0, 0, 2.0_dp) - &
        0.90128629936044729874_dp) < 3.0e-12_dp, &
        "Fourier facade exposes the pinned FortNum toroidal harmonic")
    call check_summary("canonical Fourier facade analytical oracle")

end program test_fourier_facade
