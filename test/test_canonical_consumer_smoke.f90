program test_canonical_consumer_smoke
    !! Small downstream-style client for the staged canonical API.
    !!
    !! This executable deliberately does not import fortfem_api or any
    !! implementation module.  Each section selects one canonical facade:
    !! geometry, Fourier-FEM, an open-boundary DtN operator,
    !! structure-preserving time stepping, and interchange samples.
    use check, only: check_condition, check_summary
    use fortfem_core, only: toroidal_point_to_cartesian
    use fortfem_fourier, only: &
        evaluate_fourier_mode_expansion, fourier_mode_registry_t, &
        initialize_fourier_mode_registry
    use fortfem_interop, only: &
        compare_interchange_samples, initialize_interchange_samples, &
        interchange_sample_set_t, validate_interchange_samples
    use fortfem_kinds, only: dp
    use fortfem_time, only: advance_mixed_wave_midpoint
    use fortfem_boundary, only: apply_planar_helmholtz_dtn
    use fortsparse, only: fortsparse_status_t
    implicit none

    call check_core_consumer()
    call check_fourier_consumer()
    call check_boundary_consumer()
    call check_mixed_wave_consumer()
    call check_interchange_consumer()
    call check_summary("canonical downstream consumer clients")

contains

    subroutine check_core_consumer()
        real(dp), parameter :: scale = 2.0_dp, eta = 0.65_dp
        real(dp), parameter :: theta = 0.37_dp, phi = -0.41_dp
        real(dp) :: point(3), denominator, radial, vertical, expected(3)

        denominator = cosh(eta) - cos(theta)
        radial = scale*sinh(eta)/denominator
        vertical = scale*sin(theta)/denominator
        expected = [radial*cos(phi), radial*sin(phi), vertical]
        call toroidal_point_to_cartesian(scale, eta, theta, phi, point)
        call check_condition(maxval(abs(point - expected)) < 2.0e-14_dp, &
            "downstream client imports the core toroidal map and gets its analytical value")
    end subroutine check_core_consumer

    subroutine check_fourier_consumer()
        integer, parameter :: mode_count = 1
        integer, parameter :: poloidal(mode_count) = [1]
        integer, parameter :: toroidal(mode_count) = [0]
        integer, parameter :: radial_power(mode_count) = [1]
        real(dp), parameter :: normalization(mode_count) = [1.0_dp]
        real(dp), parameter :: radius = 0.45_dp, theta = 0.31_dp
        complex(dp), parameter :: coefficient = cmplx(1.2_dp, -0.4_dp, dp)
        complex(dp) :: coefficients(mode_count), value, radial_derivative
        complex(dp) :: theta_derivative, phi_derivative, phase, expected
        type(fourier_mode_registry_t) :: registry
        type(fortsparse_status_t) :: status

        coefficients = [coefficient]
        call initialize_fourier_mode_registry( &
            registry, poloidal, toroidal, 1, 0.0_dp, 0.0_dp, .false., &
            radial_powers=radial_power, normalization=normalization, status=status)
        call check_condition(status%code == 0, &
            "downstream Fourier client initializes a canonical mode registry")
        call evaluate_fourier_mode_expansion( &
            registry, coefficients, radius, theta, 0.0_dp, value, &
            radial_derivative, theta_derivative, phi_derivative, status)
        phase = cmplx(cos(theta), sin(theta), dp)
        expected = coefficient*radius*phase
        call check_condition(status%code == 0 .and. abs(value - expected) < 2.0e-13_dp, &
            "downstream Fourier client agrees with the independent single-mode oracle")
        call check_condition(abs(radial_derivative - coefficient*phase) < 2.0e-13_dp .and. &
            abs(theta_derivative - cmplx(0.0_dp, 1.0_dp, dp)*expected) < 2.0e-13_dp .and. &
            abs(phi_derivative) < 2.0e-13_dp, &
            "Fourier client exposes the expected radial, poloidal, and toroidal derivatives")
    end subroutine check_fourier_consumer

    subroutine check_boundary_consumer()
        integer, parameter :: sample_count = 8
        real(dp), parameter :: wavenumber = 2.5_dp, period = 2.0_dp*acos(-1.0_dp)
        complex(dp) :: trace(sample_count), normal_derivative(sample_count)

        trace = cmplx(1.0_dp, 0.0_dp, dp)
        call apply_planar_helmholtz_dtn( &
            trace, wavenumber, period, normal_derivative)
        call check_condition(maxval(abs(normal_derivative - &
            cmplx(0.0_dp, wavenumber, dp))) < 2.0e-12_dp, &
            "downstream boundary client maps the constant trace to the exact DtN mode")
    end subroutine check_boundary_consumer

    subroutine check_mixed_wave_consumer()
        integer, parameter :: step_count = 20
        real(dp), parameter :: time_step = 0.01_dp, frequency = 1.7_dp
        real(dp) :: mass_q(1, 1), mass_v(1, 1), coupling(1, 1)
        real(dp) :: q(1), v(1), q_initial(1), v_initial(1)
        real(dp) :: q_exact, v_exact, energy_initial, energy
        real(dp) :: time, maximum_error, maximum_energy_drift
        integer :: step
        type(fortsparse_status_t) :: status

        mass_q = 1.0_dp
        mass_v = 1.0_dp
        coupling = frequency
        q = [1.0_dp]
        v = [0.2_dp]
        q_initial = q
        v_initial = v
        energy_initial = 0.5_dp*(q(1)**2 + v(1)**2)
        maximum_error = 0.0_dp
        maximum_energy_drift = 0.0_dp
        do step = 1, step_count
            call advance_mixed_wave_midpoint( &
                mass_q, mass_v, coupling, time_step, q, v, status)
            call check_condition(status%code == 0, &
                "downstream mixed-wave client accepts compatible mass and coupling blocks")
            time = real(step, dp)*time_step
            q_exact = q_initial(1)*cos(frequency*time) - &
                v_initial(1)*sin(frequency*time)
            v_exact = v_initial(1)*cos(frequency*time) + &
                q_initial(1)*sin(frequency*time)
            maximum_error = max(maximum_error, abs(q(1) - q_exact), &
                abs(v(1) - v_exact))
            energy = 0.5_dp*(q(1)**2 + v(1)**2)
            maximum_energy_drift = max(maximum_energy_drift, &
                abs(energy - energy_initial))
        end do
        call check_condition(maximum_error < 2.0e-4_dp, &
            "mixed-wave client follows the independent harmonic oscillator oracle")
        call check_condition(maximum_energy_drift < 2.0e-13_dp, &
            "mixed-wave client preserves its quadratic energy")

        do step = 1, step_count
            call advance_mixed_wave_midpoint( &
                mass_q, mass_v, coupling, -time_step, q, v, status)
        end do
        call check_condition(maxval(abs(q - q_initial)) < 2.0e-13_dp .and. &
            maxval(abs(v - v_initial)) < 2.0e-13_dp, &
            "mixed-wave client is exactly reversible under signed time-step reversal")
    end subroutine check_mixed_wave_consumer

    subroutine check_interchange_consumer()
        real(dp) :: coordinates(1, 3), reference_values(1, 3)
        real(dp) :: candidate_values(1, 3), weights(3)
        real(dp) :: absolute_error, relative_error, maximum_error
        type(interchange_sample_set_t) :: reference, candidate
        integer :: status

        coordinates = reshape([0.0_dp, 0.5_dp, 1.0_dp], shape(coordinates))
        reference_values = reshape([2.0_dp, 3.5_dp, 5.0_dp], &
            shape(reference_values))
        candidate_values = reference_values
        candidate_values(1, 2) = candidate_values(1, 2) + 0.25_dp
        weights = [1.0_dp, 2.0_dp, 1.0_dp]
        call initialize_interchange_samples( &
            reference, coordinates, reference_values, weights, "analytic", &
            "u(x)=2+3x", status)
        call check_condition(status == 0 .and. &
            validate_interchange_samples(reference, status), &
            "downstream interchange client validates analytical samples and provenance")
        call initialize_interchange_samples( &
            candidate, coordinates, candidate_values, weights, "consumer", &
            "perturbed analytic", status)
        call compare_interchange_samples( &
            reference, candidate, 1.0e-14_dp, absolute_error, relative_error, &
            maximum_error, status)
        call check_condition(status == 0 .and. &
            abs(absolute_error - sqrt(2.0_dp)*0.25_dp) < 2.0e-14_dp .and. &
            maximum_error == 0.25_dp .and. relative_error > 0.0_dp, &
            "interchange client reproduces the independent weighted sample error")
    end subroutine check_interchange_consumer

end program test_canonical_consumer_smoke
