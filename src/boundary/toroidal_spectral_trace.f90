module fortfem_toroidal_spectral_trace
    !! Fourier/toroidal-harmonic trace map for smooth constant-eta surfaces.
    !!
    !! This is a neutral NESTOR-like spectral building block: callers provide
    !! modal coefficients and choose the regular P branch or decaying Q branch.
    !! No equilibrium profile, coil model, or toroidal field convention is
    !! selected here.
    use fortfem_kinds, only: dp
    use fortnum_special_toroidal, only: &
        toroidal_p, toroidal_p_derivative, toroidal_q, toroidal_q_derivative
    implicit none
    private

    complex(dp), parameter :: imaginary_unit = cmplx(0.0_dp, 1.0_dp, dp)

    public :: evaluate_toroidal_spectral_trace
    public :: evaluate_toroidal_spectral_trace_jvp
    public :: evaluate_toroidal_spectral_trace_vjp
    public :: evaluate_toroidal_spectral_trace_grid
    public :: evaluate_toroidal_spectral_trace_grid_jvp
    public :: evaluate_toroidal_spectral_trace_grid_vjp

contains

    subroutine evaluate_toroidal_spectral_trace( &
            degree_indices, orders, coefficients, scale, eta, theta, phi, &
            use_second_kind, value, normal_derivative, status)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: scale, eta, theta, phi
        logical, intent(in) :: use_second_kind
        complex(dp), intent(out) :: value, normal_derivative
        integer, intent(out) :: status

        complex(dp) :: mode_value, mode_normal, mode_value_dot, mode_normal_dot
        integer :: mode

        value = cmplx(0.0_dp, 0.0_dp, dp)
        normal_derivative = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_inputs( &
            degree_indices, orders, coefficients, scale, eta)) return
        do mode = 1, size(coefficients)
            call evaluate_mode( &
                degree_indices(mode), orders(mode), scale, eta, theta, phi, &
                use_second_kind, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                mode_value, mode_normal, mode_value_dot, mode_normal_dot, status)
            if (status /= 0) return
            value = value + coefficients(mode)*mode_value
            normal_derivative = normal_derivative + coefficients(mode)*mode_normal
        end do
        status = 0
    end subroutine evaluate_toroidal_spectral_trace

    subroutine evaluate_toroidal_spectral_trace_jvp( &
            degree_indices, orders, coefficients, scale, eta, theta, phi, &
            use_second_kind, coefficients_dot, scale_dot, eta_dot, theta_dot, &
            phi_dot, value_dot, normal_derivative_dot, status)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: coefficients(:), coefficients_dot(:)
        real(dp), intent(in) :: scale, eta, theta, phi
        logical, intent(in) :: use_second_kind
        real(dp), intent(in) :: scale_dot, eta_dot, theta_dot, phi_dot
        complex(dp), intent(out) :: value_dot, normal_derivative_dot
        integer, intent(out) :: status

        complex(dp) :: mode_value, mode_normal, mode_value_dot, mode_normal_dot
        integer :: mode

        value_dot = cmplx(0.0_dp, 0.0_dp, dp)
        normal_derivative_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_inputs( &
            degree_indices, orders, coefficients, scale, eta)) return
        if (size(coefficients_dot) /= size(coefficients)) return
        do mode = 1, size(coefficients)
            call evaluate_mode( &
                degree_indices(mode), orders(mode), scale, eta, theta, phi, &
                use_second_kind, scale_dot, eta_dot, theta_dot, phi_dot, 0.0_dp, &
                mode_value, mode_normal, mode_value_dot, mode_normal_dot, status)
            if (status /= 0) return
            value_dot = value_dot + coefficients_dot(mode)*mode_value + &
                coefficients(mode)*mode_value_dot
            normal_derivative_dot = normal_derivative_dot + &
                coefficients_dot(mode)*mode_normal + &
                coefficients(mode)*mode_normal_dot
        end do
        status = 0
    end subroutine evaluate_toroidal_spectral_trace_jvp

    subroutine evaluate_toroidal_spectral_trace_grid( &
            degree_indices, orders, coefficients, scale, eta, theta, phi, &
            use_second_kind, values, normal_derivatives, status)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: scale, eta(:), theta(:), phi(:)
        logical, intent(in) :: use_second_kind
        complex(dp), intent(out) :: values(:), normal_derivatives(:)
        integer, intent(out) :: status

        integer :: point

        values = cmplx(0.0_dp, 0.0_dp, dp)
        normal_derivatives = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (size(eta) < 1 .or. size(theta) /= size(eta) .or. &
            size(phi) /= size(eta) .or. size(values) /= size(eta) .or. &
            size(normal_derivatives) /= size(eta)) return
        do point = 1, size(eta)
            call evaluate_toroidal_spectral_trace( &
                degree_indices, orders, coefficients, scale, eta(point), &
                theta(point), phi(point), use_second_kind, values(point), &
                normal_derivatives(point), status)
            if (status /= 0) return
        end do
        status = 0
    end subroutine evaluate_toroidal_spectral_trace_grid

    subroutine evaluate_toroidal_spectral_trace_grid_jvp( &
            degree_indices, orders, coefficients, scale, eta, theta, phi, &
            use_second_kind, coefficients_dot, scale_dot, eta_dot, theta_dot, &
            phi_dot, values_dot, normal_derivatives_dot, status)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: coefficients(:), coefficients_dot(:)
        real(dp), intent(in) :: scale, eta(:), theta(:), phi(:)
        logical, intent(in) :: use_second_kind
        real(dp), intent(in) :: scale_dot, eta_dot(:), theta_dot(:), phi_dot(:)
        complex(dp), intent(out) :: values_dot(:), normal_derivatives_dot(:)
        integer, intent(out) :: status

        integer :: point

        values_dot = cmplx(0.0_dp, 0.0_dp, dp)
        normal_derivatives_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (size(eta) < 1 .or. size(theta) /= size(eta) .or. &
            size(phi) /= size(eta) .or. size(eta_dot) /= size(eta) .or. &
            size(theta_dot) /= size(eta) .or. size(phi_dot) /= size(eta) .or. &
            size(values_dot) /= size(eta) .or. &
            size(normal_derivatives_dot) /= size(eta)) return
        do point = 1, size(eta)
            call evaluate_toroidal_spectral_trace_jvp( &
                degree_indices, orders, coefficients, scale, eta(point), &
                theta(point), phi(point), use_second_kind, coefficients_dot, &
                scale_dot, eta_dot(point), theta_dot(point), phi_dot(point), &
                values_dot(point), normal_derivatives_dot(point), status)
            if (status /= 0) return
        end do
        status = 0
    end subroutine evaluate_toroidal_spectral_trace_grid_jvp

    subroutine evaluate_toroidal_spectral_trace_grid_vjp( &
            degree_indices, orders, coefficients, scale, eta, theta, phi, &
            use_second_kind, values_bar, normal_derivatives_bar, coefficients_bar, &
            scale_bar, eta_bar, theta_bar, phi_bar, status)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: coefficients(:), values_bar(:)
        complex(dp), intent(in) :: normal_derivatives_bar(:)
        real(dp), intent(in) :: scale, eta(:), theta(:), phi(:)
        logical, intent(in) :: use_second_kind
        complex(dp), allocatable, intent(out) :: coefficients_bar(:)
        real(dp), intent(out) :: scale_bar, eta_bar(:), theta_bar(:), phi_bar(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: local_coefficients_bar(:)
        real(dp) :: local_scale_bar, local_eta_bar, local_theta_bar, local_phi_bar
        integer :: point

        scale_bar = 0.0_dp
        eta_bar = 0.0_dp
        theta_bar = 0.0_dp
        phi_bar = 0.0_dp
        if (allocated(coefficients_bar)) deallocate(coefficients_bar)
        status = 1
        if (size(eta) < 1 .or. size(theta) /= size(eta) .or. &
            size(phi) /= size(eta) .or. size(values_bar) /= size(eta) .or. &
            size(normal_derivatives_bar) /= size(eta) .or. &
            size(eta_bar) /= size(eta) .or. size(theta_bar) /= size(eta) .or. &
            size(phi_bar) /= size(eta)) return
        allocate(coefficients_bar(size(coefficients)))
        coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        do point = 1, size(eta)
            call evaluate_toroidal_spectral_trace_vjp( &
                degree_indices, orders, coefficients, scale, eta(point), &
                theta(point), phi(point), use_second_kind, values_bar(point), &
                normal_derivatives_bar(point), local_coefficients_bar, &
                local_scale_bar, local_eta_bar, local_theta_bar, local_phi_bar, status)
            if (status /= 0) return
            coefficients_bar = coefficients_bar + local_coefficients_bar
            scale_bar = scale_bar + local_scale_bar
            eta_bar(point) = local_eta_bar
            theta_bar(point) = local_theta_bar
            phi_bar(point) = local_phi_bar
        end do
        status = 0
    end subroutine evaluate_toroidal_spectral_trace_grid_vjp

    subroutine evaluate_toroidal_spectral_trace_vjp( &
            degree_indices, orders, coefficients, scale, eta, theta, phi, &
            use_second_kind, value_bar, normal_derivative_bar, coefficients_bar, &
            scale_bar, eta_bar, theta_bar, phi_bar, status)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: coefficients(:), value_bar, normal_derivative_bar
        real(dp), intent(in) :: scale, eta, theta, phi
        logical, intent(in) :: use_second_kind
        complex(dp), allocatable, intent(out) :: coefficients_bar(:)
        real(dp), intent(out) :: scale_bar, eta_bar, theta_bar, phi_bar
        integer, intent(out) :: status

        complex(dp) :: mode_value, mode_normal, mode_value_dot, mode_normal_dot
        integer :: mode

        scale_bar = 0.0_dp
        eta_bar = 0.0_dp
        theta_bar = 0.0_dp
        phi_bar = 0.0_dp
        if (allocated(coefficients_bar)) deallocate(coefficients_bar)
        status = 1
        if (.not. valid_inputs( &
            degree_indices, orders, coefficients, scale, eta)) return
        allocate(coefficients_bar(size(coefficients)))
        coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        do mode = 1, size(coefficients)
            call evaluate_mode( &
                degree_indices(mode), orders(mode), scale, eta, theta, phi, &
                use_second_kind, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                mode_value, mode_normal, mode_value_dot, mode_normal_dot, status)
            if (status /= 0) return
            coefficients_bar(mode) = value_bar*conjg(mode_value) + &
                normal_derivative_bar*conjg(mode_normal)
            call mode_direction( &
                degree_indices(mode), orders(mode), scale, eta, theta, phi, &
                use_second_kind, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                mode_value_dot, mode_normal_dot, status)
            if (status /= 0) return
            scale_bar = scale_bar + real(conjg(value_bar)*coefficients(mode)* &
                mode_value_dot + conjg(normal_derivative_bar)*coefficients(mode)* &
                mode_normal_dot, dp)
            call mode_direction( &
                degree_indices(mode), orders(mode), scale, eta, theta, phi, &
                use_second_kind, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                mode_value_dot, mode_normal_dot, status)
            if (status /= 0) return
            eta_bar = eta_bar + real(conjg(value_bar)*coefficients(mode)* &
                mode_value_dot + conjg(normal_derivative_bar)*coefficients(mode)* &
                mode_normal_dot, dp)
            call mode_direction( &
                degree_indices(mode), orders(mode), scale, eta, theta, phi, &
                use_second_kind, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                mode_value_dot, mode_normal_dot, status)
            if (status /= 0) return
            theta_bar = theta_bar + real(conjg(value_bar)*coefficients(mode)* &
                mode_value_dot + conjg(normal_derivative_bar)*coefficients(mode)* &
                mode_normal_dot, dp)
            call mode_direction( &
                degree_indices(mode), orders(mode), scale, eta, theta, phi, &
                use_second_kind, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                mode_value_dot, mode_normal_dot, status)
            if (status /= 0) return
            phi_bar = phi_bar + real(conjg(value_bar)*coefficients(mode)* &
                mode_value_dot + conjg(normal_derivative_bar)*coefficients(mode)* &
                mode_normal_dot, dp)
        end do
        status = 0
    end subroutine evaluate_toroidal_spectral_trace_vjp

    subroutine mode_direction( &
            degree_index, order, scale, eta, theta, phi, use_second_kind, &
            scale_dot, eta_dot, theta_dot, phi_dot, unused, value_dot, &
            normal_dot, status)
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: scale, eta, theta, phi
        logical, intent(in) :: use_second_kind
        real(dp), intent(in) :: scale_dot, eta_dot, theta_dot, phi_dot, unused
        complex(dp), intent(out) :: value_dot, normal_dot
        integer, intent(out) :: status
        complex(dp) :: value, normal

        call evaluate_mode( &
            degree_index, order, scale, eta, theta, phi, use_second_kind, &
            scale_dot, eta_dot, theta_dot, phi_dot, unused, value, normal, &
            value_dot, normal_dot, status)
    end subroutine mode_direction

    subroutine evaluate_mode( &
            degree_index, order, scale, eta, theta, phi, use_second_kind, &
            scale_dot, eta_dot, theta_dot, phi_dot, unused, value, normal, &
            value_dot, normal_dot, status)
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: scale, eta, theta, phi
        logical, intent(in) :: use_second_kind
        real(dp), intent(in) :: scale_dot, eta_dot, theta_dot, phi_dot, unused
        complex(dp), intent(out) :: value, normal, value_dot, normal_dot
        integer, intent(out) :: status

        real(dp) :: denominator, denominator_dot, degree, radial, radial_dot
        real(dp) :: radial_first, radial_second, radial_first_dot
        real(dp) :: square_root, square_root_dot, x, x_dot, sine_eta
        real(dp) :: coefficient, coefficient_dot
        real(dp) :: eta_factor, eta_factor_dot, phase_angle_dot
        real(dp) :: normal_factor, normal_factor_dot
        complex(dp) :: phase, phase_dot

        value = cmplx(0.0_dp, 0.0_dp, dp)
        normal = cmplx(0.0_dp, 0.0_dp, dp)
        value_dot = cmplx(0.0_dp, 0.0_dp, dp)
        normal_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (degree_index < 0 .or. order < 0 .or. scale <= 0.0_dp .or. &
            eta <= 0.0_dp) return
        x = cosh(eta)
        sine_eta = sinh(eta)
        denominator = x - cos(theta)
        if (denominator <= tiny(1.0_dp)) return
        square_root = sqrt(denominator)
        degree = real(degree_index, dp) - 0.5_dp
        if (use_second_kind) then
            radial = toroidal_q(degree_index, order, x)
            radial_first = toroidal_q_derivative(degree_index, order, x)
        else
            radial = toroidal_p(degree_index, order, x)
            radial_first = toroidal_p_derivative(degree_index, order, x)
        end if
        if (radial /= radial .or. radial_first /= radial_first) return
        radial_second = ((degree*(degree + 1.0_dp) + &
            real(order*order, dp)/(x*x - 1.0_dp))*radial - &
            2.0_dp*x*radial_first)/(x*x - 1.0_dp)
        phase = exp(imaginary_unit*real( &
            degree_index, dp)*theta + imaginary_unit*real(order, dp)*phi)
        coefficient = radial/(2.0_dp*square_root) + square_root*radial_first
        eta_factor = sine_eta*coefficient
        value = square_root*radial*phase
        normal_factor = -denominator/scale
        normal = normal_factor*eta_factor*phase

        denominator_dot = sine_eta*eta_dot + sin(theta)*theta_dot
        square_root_dot = denominator_dot/(2.0_dp*square_root)
        x_dot = sine_eta*eta_dot
        radial_dot = radial_first*x_dot
        radial_first_dot = radial_second*x_dot
        phase_angle_dot = real(degree_index, dp)*theta_dot + &
            real(order, dp)*phi_dot
        phase_dot = imaginary_unit*phase_angle_dot*phase
        coefficient_dot = radial_dot/(2.0_dp*square_root) - &
            radial*square_root_dot/(2.0_dp*denominator) + &
            square_root_dot*radial_first + square_root*radial_first_dot
        eta_factor_dot = cosh(eta)*eta_dot*coefficient + &
            sine_eta*coefficient_dot
        value_dot = (square_root_dot*radial + square_root*radial_dot)*phase + &
            square_root*radial*phase_dot
        normal_factor_dot = -(denominator_dot/scale - &
            denominator*scale_dot/scale**2)
        normal_dot = normal_factor_dot*eta_factor*phase + &
            normal_factor*eta_factor_dot*phase + &
            normal_factor*eta_factor*phase_dot
        status = 0
    end subroutine evaluate_mode

    logical function valid_inputs( &
            degree_indices, orders, coefficients, scale, eta) result(valid)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: coefficients(:)
        real(dp), intent(in) :: scale, eta

        valid = size(degree_indices) > 0 .and. &
            size(orders) == size(degree_indices) .and. &
            size(coefficients) == size(degree_indices) .and. &
            all(degree_indices >= 0) .and. all(orders >= 0) .and. &
            scale > 0.0_dp .and. eta > 0.0_dp
    end function valid_inputs

end module fortfem_toroidal_spectral_trace
