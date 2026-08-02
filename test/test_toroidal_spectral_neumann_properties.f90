program test_toroidal_spectral_neumann_properties
    use, intrinsic :: iso_fortran_env, only: int32
    use check, only: check_condition, check_property, check_summary, &
        property_random_integer, property_random_unit, property_rng_initialize, &
        property_rng_t
    use fortfem_api, only: &
        solve_toroidal_spectral_neumann, solve_toroidal_spectral_neumann_jvp, &
        solve_toroidal_spectral_neumann_vjp, toroidal_p, toroidal_p_derivative, &
        toroidal_q, toroidal_q_derivative
    use fortfem_kinds, only: dp
    implicit none

    logical :: all_passed
    integer(int32) :: first_failed_seed, shrunk_seed

    call check_property( &
        "seeded toroidal Neumann modal division and derivative identities", &
        864209_int32, 12, neumann_case, all_passed, first_failed_seed, &
        shrunk_seed)
    call check_condition(all_passed .and. first_failed_seed == 0_int32 .and. &
        shrunk_seed == 0_int32, &
        "random toroidal Neumann property reports no failure seed")
    call check_summary("random toroidal Neumann properties")
    if (.not. all_passed) error stop 1

contains

    logical function neumann_case(case_seed)
        integer(int32), intent(in) :: case_seed
        real(dp), parameter :: epsilon_fd = 1.0e-5_dp
        real(dp), parameter :: resonance_tolerance = 1.0e-12_dp
        integer, parameter :: maximum_modes = 4
        type(property_rng_t) :: rng
        integer, allocatable :: degree_indices(:), orders(:)
        integer, allocatable :: zero_degrees(:), zero_orders(:)
        complex(dp), allocatable :: normal_coefficients(:)
        complex(dp), allocatable :: normal_coefficients_dot(:)
        complex(dp), allocatable :: potential(:), potential_dot(:)
        complex(dp), allocatable :: potential_plus(:), potential_minus(:)
        complex(dp), allocatable :: potential_bar(:), normal_coefficients_bar(:)
        real(dp) :: scale, eta, scale_dot, eta_dot
        real(dp) :: scale_bar, eta_bar, lhs, rhs
        real(dp) :: derivative_error, adjoint_error
        integer :: mode_count, mode, status, oracle_status
        logical :: use_second_kind

        neumann_case = .false.
        call property_rng_initialize(rng, case_seed)
        mode_count = property_random_integer(rng, 1, maximum_modes)
        allocate(degree_indices(mode_count), orders(mode_count))
        allocate(zero_degrees(mode_count), zero_orders(mode_count))
        allocate(normal_coefficients(mode_count), normal_coefficients_dot(mode_count))
        allocate(potential(mode_count), potential_dot(mode_count))
        allocate(potential_plus(mode_count), potential_minus(mode_count))
        allocate(potential_bar(mode_count), normal_coefficients_bar(mode_count))

        do mode = 1, mode_count
            ! FortNum's Q derivative is defined for positive degree and
            ! order <= degree; the same valid set is used for both branches.
            degree_indices(mode) = property_random_integer(rng, 1, 3)
            orders(mode) = property_random_integer(rng, 0, degree_indices(mode))
            normal_coefficients(mode) = random_complex(rng)
            normal_coefficients_dot(mode) = 0.25_dp*random_complex(rng)
            potential_bar(mode) = random_complex(rng)
        end do
        scale = 1.1_dp + 0.9_dp*property_random_unit(rng)
        eta = 0.6_dp + 0.9_dp*property_random_unit(rng)
        scale_dot = 0.15_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
        eta_dot = 0.15_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
        use_second_kind = property_random_integer(rng, 0, 1) == 1

        call modal_division_oracle( &
            degree_indices, orders, normal_coefficients, scale, eta, &
            use_second_kind, potential, oracle_status)
        if (oracle_status /= 0) return
        call solve_toroidal_spectral_neumann( &
            degree_indices, orders, normal_coefficients, scale, eta, &
            use_second_kind, .true., resonance_tolerance, potential_dot, status)
        if (status /= 0) return
        if (maxval(abs(potential_dot - potential)) > 3.0e-12_dp) return

        call solve_toroidal_spectral_neumann_jvp( &
            degree_indices, orders, normal_coefficients, scale, eta, &
            use_second_kind, .true., resonance_tolerance, normal_coefficients_dot, &
            scale_dot, eta_dot, potential_dot, status)
        if (status /= 0) return
        call solve_toroidal_spectral_neumann( &
            degree_indices, orders, normal_coefficients + epsilon_fd* &
            normal_coefficients_dot, scale + epsilon_fd*scale_dot, &
            eta + epsilon_fd*eta_dot, use_second_kind, .true., &
            resonance_tolerance, potential_plus, status)
        if (status /= 0) return
        call solve_toroidal_spectral_neumann( &
            degree_indices, orders, normal_coefficients - epsilon_fd* &
            normal_coefficients_dot, scale - epsilon_fd*scale_dot, &
            eta - epsilon_fd*eta_dot, use_second_kind, .true., &
            resonance_tolerance, potential_minus, status)
        if (status /= 0) return
        derivative_error = maxval(abs(potential_dot - &
            (potential_plus - potential_minus)/(2.0_dp*epsilon_fd)))
        if (derivative_error > 3.0e-7_dp) return

        call solve_toroidal_spectral_neumann_vjp( &
            degree_indices, orders, normal_coefficients, scale, eta, &
            use_second_kind, .true., resonance_tolerance, potential_bar, &
            normal_coefficients_bar, scale_bar, eta_bar, status)
        if (status /= 0) return
        lhs = real(sum(conjg(potential_bar)*potential_dot), dp)
        rhs = real(sum(conjg(normal_coefficients_bar)*normal_coefficients_dot), dp) + &
            scale_bar*scale_dot + eta_bar*eta_dot
        adjoint_error = abs(lhs - rhs)
        if (adjoint_error > 5.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs))) return

        zero_degrees = degree_indices
        zero_orders = orders
        zero_degrees(1) = 0
        zero_orders(1) = 0
        call solve_toroidal_spectral_neumann( &
            zero_degrees, zero_orders, normal_coefficients, scale, eta, .false., &
            .false., resonance_tolerance, potential, status)
        if (status == 0) return

        call solve_toroidal_spectral_neumann( &
            degree_indices, orders, normal_coefficients, scale, eta, &
            use_second_kind, .true., 1.0e12_dp, potential, status)
        if (status == 0) return

        neumann_case = .true.
    end function neumann_case

    function random_complex(rng) result(value)
        type(property_rng_t), intent(inout) :: rng
        complex(dp) :: value

        value = cmplx(2.0_dp*property_random_unit(rng) - 1.0_dp, &
            2.0_dp*property_random_unit(rng) - 1.0_dp, dp)
    end function random_complex

    subroutine modal_division_oracle( &
            degree_indices, orders, normal_coefficients, scale, eta, &
            use_second_kind, potential, status)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: normal_coefficients(:)
        real(dp), intent(in) :: scale, eta
        logical, intent(in) :: use_second_kind
        complex(dp), intent(out) :: potential(:)
        integer, intent(out) :: status
        real(dp) :: denominator, square_root, radial, radial_first
        real(dp) :: mode_normal, x
        integer :: mode

        potential = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (size(degree_indices) < 1 .or. size(orders) /= size(degree_indices) .or. &
            size(normal_coefficients) /= size(degree_indices) .or. &
            size(potential) /= size(normal_coefficients) .or. &
            any(degree_indices < 1) .or. any(orders < 0) .or. &
            any(orders > degree_indices) .or. scale <= 0.0_dp .or. eta <= 0.0_dp) return
        x = cosh(eta)
        denominator = x - 1.0_dp
        square_root = sqrt(denominator)
        do mode = 1, size(normal_coefficients)
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
            mode_normal = -denominator/scale*sinh(eta)* &
                (radial/(2.0_dp*square_root) + square_root*radial_first)
            if (abs(mode_normal) <= tiny(1.0_dp)) return
            potential(mode) = normal_coefficients(mode)/mode_normal
        end do
        status = 0
    end subroutine modal_division_oracle

end program test_toroidal_spectral_neumann_properties
