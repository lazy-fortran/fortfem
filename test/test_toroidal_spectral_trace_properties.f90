program test_toroidal_spectral_trace_properties
    use, intrinsic :: iso_fortran_env, only: int32
    use check, only: check_condition, check_property, check_summary, &
        property_random_integer, property_random_unit, property_rng_initialize, &
        property_rng_t
    use fortfem_api, only: &
        evaluate_toroidal_spectral_trace, &
        evaluate_toroidal_spectral_trace_grid, &
        evaluate_toroidal_spectral_trace_grid_jvp, &
        evaluate_toroidal_spectral_trace_jvp, toroidal_p, toroidal_p_derivative, &
        toroidal_q, toroidal_q_derivative
    use fortfem_kinds, only: dp
    implicit none

    logical :: all_passed
    integer(int32) :: first_failed_seed, shrunk_seed

    call check_property( &
        "seeded toroidal P/Q spectral trace modal and derivative identities", &
        97531_int32, 12, spectral_case, all_passed, first_failed_seed, &
        shrunk_seed)
    call check_condition(all_passed .and. first_failed_seed == 0_int32 .and. &
        shrunk_seed == 0_int32, &
        "random toroidal spectral property reports no failure seed")
    call check_summary("random toroidal spectral trace properties")
    if (.not. all_passed) error stop 1

contains

    logical function spectral_case(case_seed)
        integer(int32), intent(in) :: case_seed
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp), parameter :: epsilon_fd = 1.0e-5_dp
        integer, parameter :: maximum_modes = 4
        integer, parameter :: maximum_points = 4
        type(property_rng_t) :: rng
        integer, allocatable :: degree_indices(:), orders(:), bad_degrees(:)
        complex(dp), allocatable :: coefficients(:), coefficients_dot(:)
        complex(dp), allocatable :: values(:), normal_derivatives(:)
        complex(dp), allocatable :: values_dot(:), normal_derivatives_dot(:)
        complex(dp), allocatable :: values_plus(:), values_minus(:)
        complex(dp), allocatable :: normals_plus(:), normals_minus(:)
        real(dp), allocatable :: eta(:), theta(:), phi(:)
        real(dp), allocatable :: eta_dot(:), theta_dot(:), phi_dot(:)
        real(dp) :: scale, scale_dot
        real(dp) :: value_error, normal_error
        real(dp) :: value_single_error, normal_single_error
        complex(dp) :: value_single, normal_single
        complex(dp) :: value_plus, value_minus, normal_plus, normal_minus
        complex(dp) :: expected_value, expected_normal
        integer :: mode_count, point_count, mode, point, status
        integer :: oracle_status
        logical :: use_second_kind

        spectral_case = .false.
        call property_rng_initialize(rng, case_seed)
        mode_count = property_random_integer(rng, 1, maximum_modes)
        point_count = property_random_integer(rng, 2, maximum_points)
        allocate(degree_indices(mode_count), orders(mode_count), &
            coefficients(mode_count), coefficients_dot(mode_count))
        allocate(eta(point_count), theta(point_count), phi(point_count))
        allocate(eta_dot(point_count), theta_dot(point_count), phi_dot(point_count))
        allocate(values(point_count), normal_derivatives(point_count))
        allocate(values_dot(point_count), normal_derivatives_dot(point_count))
        allocate(values_plus(point_count), values_minus(point_count))
        allocate(normals_plus(point_count), normals_minus(point_count))

        do mode = 1, mode_count
            degree_indices(mode) = property_random_integer(rng, 1, 3)
            orders(mode) = property_random_integer(rng, 0, degree_indices(mode))
            coefficients(mode) = random_complex(rng)
            coefficients_dot(mode) = 0.2_dp*random_complex(rng)
        end do
        do point = 1, point_count
            eta(point) = 0.45_dp + 1.0_dp*property_random_unit(rng)
            theta(point) = pi*(2.0_dp*property_random_unit(rng) - 1.0_dp)
            phi(point) = pi*(2.0_dp*property_random_unit(rng) - 1.0_dp)
            eta_dot(point) = 0.12_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
            theta_dot(point) = 0.12_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
            phi_dot(point) = 0.12_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
        end do
        scale = 1.1_dp + 0.9_dp*property_random_unit(rng)
        scale_dot = 0.12_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
        use_second_kind = property_random_integer(rng, 0, 1) == 1

        call evaluate_toroidal_spectral_trace_grid( &
            degree_indices, orders, coefficients, scale, eta, theta, phi, &
            use_second_kind, values, normal_derivatives, status)
        if (status /= 0) then
            return
        end if
        do point = 1, point_count
            call independent_modal_sum( &
                degree_indices, orders, coefficients, scale, eta(point), &
                theta(point), phi(point), use_second_kind, expected_value, &
                expected_normal, oracle_status)
            if (oracle_status /= 0) return
            if (abs(values(point) - expected_value) > 3.0e-12_dp .or. &
                abs(normal_derivatives(point) - expected_normal) > 3.0e-12_dp) return
        end do

        call evaluate_toroidal_spectral_trace( &
            degree_indices, orders, coefficients, scale, eta(1), theta(1), &
            phi(1), use_second_kind, value_single, normal_single, status)
        if (status /= 0) return
        value_single_error = abs(value_single - values(1))
        normal_single_error = abs(normal_single - normal_derivatives(1))
        if (max(value_single_error, normal_single_error) > 3.0e-12_dp) return

        call evaluate_toroidal_spectral_trace_grid_jvp( &
            degree_indices, orders, coefficients, scale, eta, theta, phi, &
            use_second_kind, coefficients_dot, scale_dot, eta_dot, theta_dot, &
            phi_dot, values_dot, normal_derivatives_dot, status)
        if (status /= 0) return
        call evaluate_toroidal_spectral_trace_grid( &
            degree_indices, orders, coefficients + epsilon_fd*coefficients_dot, &
            scale + epsilon_fd*scale_dot, eta + epsilon_fd*eta_dot, &
            theta + epsilon_fd*theta_dot, phi + epsilon_fd*phi_dot, &
            use_second_kind, values_plus, normals_plus, status)
        if (status /= 0) return
        call evaluate_toroidal_spectral_trace_grid( &
            degree_indices, orders, coefficients - epsilon_fd*coefficients_dot, &
            scale - epsilon_fd*scale_dot, eta - epsilon_fd*eta_dot, &
            theta - epsilon_fd*theta_dot, phi - epsilon_fd*phi_dot, &
            use_second_kind, values_minus, normals_minus, status)
        if (status /= 0) return
        value_error = maxval(abs(values_dot - &
            (values_plus - values_minus)/(2.0_dp*epsilon_fd)))
        normal_error = maxval(abs(normal_derivatives_dot - &
            (normals_plus - normals_minus)/(2.0_dp*epsilon_fd)))
        if (max(value_error, normal_error) > 3.0e-7_dp) return

        call evaluate_toroidal_spectral_trace_jvp( &
            degree_indices, orders, coefficients, scale, eta(1), theta(1), &
            phi(1), use_second_kind, coefficients_dot, scale_dot, eta_dot(1), &
            theta_dot(1), phi_dot(1), value_plus, normal_plus, status)
        if (status /= 0) return
        call evaluate_toroidal_spectral_trace( &
            degree_indices, orders, coefficients + epsilon_fd*coefficients_dot, &
            scale + epsilon_fd*scale_dot, eta(1) + epsilon_fd*eta_dot(1), &
            theta(1) + epsilon_fd*theta_dot(1), phi(1) + epsilon_fd*phi_dot(1), &
            use_second_kind, value_single, normal_single, status)
        if (status /= 0) return
        call evaluate_toroidal_spectral_trace( &
            degree_indices, orders, coefficients - epsilon_fd*coefficients_dot, &
            scale - epsilon_fd*scale_dot, eta(1) - epsilon_fd*eta_dot(1), &
            theta(1) - epsilon_fd*theta_dot(1), phi(1) - epsilon_fd*phi_dot(1), &
            use_second_kind, value_minus, normal_minus, status)
        if (status /= 0) return
        if (abs(value_plus - (value_single - value_minus)/ &
            (2.0_dp*epsilon_fd)) > 3.0e-7_dp .or. &
            abs(normal_plus - (normal_single - normal_minus)/ &
            (2.0_dp*epsilon_fd)) > 3.0e-7_dp) return

        bad_degrees = degree_indices
        bad_degrees(1) = -1
        call evaluate_toroidal_spectral_trace( &
            bad_degrees, orders, coefficients, scale, eta(1), theta(1), phi(1), &
            use_second_kind, value_single, normal_single, status)
        if (status == 0) return
        call evaluate_toroidal_spectral_trace_grid( &
            bad_degrees, orders, coefficients, scale, eta, theta, phi, &
            use_second_kind, values, normal_derivatives, status)
        if (status == 0) return
        call evaluate_toroidal_spectral_trace( &
            degree_indices, orders, coefficients, 0.0_dp, eta(1), theta(1), &
            phi(1), use_second_kind, value_single, normal_single, status)
        if (status == 0) return
        call evaluate_toroidal_spectral_trace( &
            degree_indices, orders, coefficients, scale, 0.0_dp, theta(1), &
            phi(1), use_second_kind, value_single, normal_single, status)
        if (status == 0) return

        spectral_case = .true.
    end function spectral_case

    function random_complex(rng) result(value)
        type(property_rng_t), intent(inout) :: rng
        complex(dp) :: value

        value = cmplx(2.0_dp*property_random_unit(rng) - 1.0_dp, &
            2.0_dp*property_random_unit(rng) - 1.0_dp, dp)
    end function random_complex

    subroutine independent_modal_sum( &
            degree_indices, orders, coefficients, scale, eta, theta, phi, &
            use_second_kind, value, normal, status)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: scale, eta, theta, phi
        logical, intent(in) :: use_second_kind
        complex(dp), intent(out) :: value, normal
        integer, intent(out) :: status
        real(dp) :: denominator, square_root, x, radial, radial_first
        complex(dp) :: phase
        integer :: mode

        value = cmplx(0.0_dp, 0.0_dp, dp)
        normal = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (size(degree_indices) < 1 .or. size(orders) /= size(degree_indices) .or. &
            size(coefficients) /= size(degree_indices) .or. &
            any(degree_indices < 0) .or. any(orders < 0) .or. scale <= 0.0_dp .or. &
            eta <= 0.0_dp) return
        x = cosh(eta)
        denominator = x - cos(theta)
        if (denominator <= tiny(1.0_dp)) return
        square_root = sqrt(denominator)
        do mode = 1, size(coefficients)
            if (use_second_kind) then
                radial = toroidal_q(degree_indices(mode), orders(mode), x)
                radial_first = toroidal_q_derivative( &
                    degree_indices(mode), orders(mode), x)
            else
                radial = toroidal_p(degree_indices(mode), orders(mode), x)
                radial_first = toroidal_p_derivative( &
                    degree_indices(mode), orders(mode), x)
            end if
            if (radial /= radial .or. radial_first /= radial_first) return
            phase = exp(cmplx(0.0_dp, real(degree_indices(mode), dp)*theta + &
                real(orders(mode), dp)*phi, dp))
            value = value + coefficients(mode)*square_root*radial*phase
            normal = normal - coefficients(mode)*denominator/scale*sinh(eta)* &
                (radial/(2.0_dp*square_root) + square_root*radial_first)*phase
        end do
        status = 0
    end subroutine independent_modal_sum

end program test_toroidal_spectral_trace_properties
