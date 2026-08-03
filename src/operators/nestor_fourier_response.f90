module fortfem_nestor_fourier_response
    !! Neutral modal response primitive for NESTOR/BIEST-like adapters.
    !!
    !! The map acts on caller-owned toroidal Fourier coefficients and returns
    !! the weighted normal radial response on a fixed eta surface,
    !!
    !!   r_mn = -(sinh(eta)/scale) R'_mn(cosh(eta)) c_mn,
    !!
    !! where R is the FortNum P or Q toroidal radial branch.  The omitted
    !! metric factor `(cosh(eta)-cos(theta))` is intentional: it is a modal
    !! radial response, while angular quadrature and Green-kernel conventions
    !! remain external to FortFEM.  This small diagonal map is therefore a
    !! reusable NESTOR/BIEST-like boundary block, not a vacuum or equilibrium
    !! solver.  All actions are fixed-topology and use the real complex-adjoint
    !! convention used by the other FortFEM response maps.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_toroidal_harmonics, only: &
        toroidal_p_derivative, toroidal_q_derivative, &
        toroidal_p_second_derivative, toroidal_q_second_derivative
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: apply_nestor_fourier_response_map
    public :: apply_nestor_fourier_response_map_jvp
    public :: apply_nestor_fourier_response_map_vjp
    public :: evaluate_nestor_fourier_reciprocity_defect

contains

    subroutine apply_nestor_fourier_response_map( &
            degree_indices, orders, source_coefficients, scale, eta, &
            use_second_kind, response_coefficients, status)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: source_coefficients(:)
        real(dp), intent(in) :: scale, eta
        logical, intent(in) :: use_second_kind
        complex(dp), intent(out) :: response_coefficients(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: factor, factor_scale, factor_eta
        integer :: mode

        response_coefficients = cmplx(0.0_dp, 0.0_dp, dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "NESTOR Fourier response received incompatible arrays")
        if (.not. valid_inputs(degree_indices, orders, source_coefficients, &
            scale, eta, response_coefficients)) return
        do mode = 1, size(source_coefficients)
            call evaluate_factor(degree_indices(mode), orders(mode), scale, eta, &
                use_second_kind, factor, factor_scale, factor_eta, status)
            if (status%code /= FORTSPARSE_OK) return
            response_coefficients(mode) = factor*source_coefficients(mode)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_nestor_fourier_response_map

    subroutine apply_nestor_fourier_response_map_jvp( &
            degree_indices, orders, source_coefficients, scale, eta, &
            use_second_kind, source_coefficients_dot, scale_dot, eta_dot, &
            response_coefficients_dot, status)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: source_coefficients(:), source_coefficients_dot(:)
        real(dp), intent(in) :: scale, eta, scale_dot, eta_dot
        logical, intent(in) :: use_second_kind
        complex(dp), intent(out) :: response_coefficients_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: factor, factor_scale, factor_eta, factor_dot
        integer :: mode

        response_coefficients_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "NESTOR Fourier response JVP received incompatible arrays")
        if (.not. valid_inputs(degree_indices, orders, source_coefficients, &
            scale, eta, response_coefficients_dot)) return
        if (size(source_coefficients_dot) /= size(source_coefficients) .or. &
            .not. finite_complex(source_coefficients_dot) .or. &
            .not. ieee_is_finite(scale_dot) .or. .not. ieee_is_finite(eta_dot)) return
        do mode = 1, size(source_coefficients)
            call evaluate_factor(degree_indices(mode), orders(mode), scale, eta, &
                use_second_kind, factor, factor_scale, factor_eta, status)
            if (status%code /= FORTSPARSE_OK) return
            factor_dot = factor_scale*scale_dot + factor_eta*eta_dot
            response_coefficients_dot(mode) = factor*source_coefficients_dot(mode) + &
                factor_dot*source_coefficients(mode)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_nestor_fourier_response_map_jvp

    subroutine apply_nestor_fourier_response_map_vjp( &
            degree_indices, orders, source_coefficients, scale, eta, &
            use_second_kind, response_coefficients_bar, source_coefficients_bar, &
            scale_bar, eta_bar, status)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: source_coefficients(:)
        real(dp), intent(in) :: scale, eta
        logical, intent(in) :: use_second_kind
        complex(dp), intent(in) :: response_coefficients_bar(:)
        complex(dp), intent(out) :: source_coefficients_bar(:)
        real(dp), intent(out) :: scale_bar, eta_bar
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: factor, factor_scale, factor_eta
        integer :: mode

        source_coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        scale_bar = 0.0_dp
        eta_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "NESTOR Fourier response VJP received incompatible arrays")
        if (.not. valid_inputs(degree_indices, orders, source_coefficients, &
            scale, eta, source_coefficients_bar)) return
        if (size(response_coefficients_bar) /= size(source_coefficients) .or. &
            .not. finite_complex(response_coefficients_bar)) return
        do mode = 1, size(source_coefficients)
            call evaluate_factor(degree_indices(mode), orders(mode), scale, eta, &
                use_second_kind, factor, factor_scale, factor_eta, status)
            if (status%code /= FORTSPARSE_OK) return
            source_coefficients_bar(mode) = conjg(cmplx(factor, 0.0_dp, dp))* &
                response_coefficients_bar(mode)
            scale_bar = scale_bar + real(conjg(response_coefficients_bar(mode))* &
                source_coefficients(mode)*factor_scale, dp)
            eta_bar = eta_bar + real(conjg(response_coefficients_bar(mode))* &
                source_coefficients(mode)*factor_eta, dp)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_nestor_fourier_response_map_vjp

    subroutine evaluate_nestor_fourier_reciprocity_defect( &
            degree_indices, orders, source_one, source_two, scale, eta, &
            use_second_kind, reciprocity_defect, status)
        !! Check the real work reciprocity of the diagonal modal response.
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: source_one(:), source_two(:)
        real(dp), intent(in) :: scale, eta
        logical, intent(in) :: use_second_kind
        real(dp), intent(out) :: reciprocity_defect
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: response_one(:), response_two(:)
        complex(dp) :: first_pairing, second_pairing
        real(dp) :: denominator

        reciprocity_defect = huge(1.0_dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "NESTOR Fourier reciprocity received incompatible arrays")
        if (size(source_one) < 1 .or. size(source_two) /= size(source_one)) return
        if (.not. finite_complex(source_one) .or. &
            .not. finite_complex(source_two)) return
        allocate(response_one(size(source_one)), response_two(size(source_one)))
        call apply_nestor_fourier_response_map( &
            degree_indices, orders, source_one, scale, eta, use_second_kind, &
            response_one, status)
        if (status%code /= FORTSPARSE_OK) return
        call apply_nestor_fourier_response_map( &
            degree_indices, orders, source_two, scale, eta, use_second_kind, &
            response_two, status)
        if (status%code /= FORTSPARSE_OK) return
        first_pairing = sum(conjg(source_one)*response_two)
        second_pairing = sum(conjg(source_two)*response_one)
        denominator = max(1.0_dp, abs(real(first_pairing, dp)), &
            abs(real(second_pairing, dp)))
        reciprocity_defect = abs(real(first_pairing-second_pairing, dp))/denominator
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_nestor_fourier_reciprocity_defect

    subroutine evaluate_factor( &
            degree_index, order, scale, eta, use_second_kind, factor, &
            factor_scale, factor_eta, status)
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: scale, eta
        logical, intent(in) :: use_second_kind
        real(dp), intent(out) :: factor, factor_scale, factor_eta
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: x, sine_eta, radial_first, radial_second

        factor = 0.0_dp
        factor_scale = 0.0_dp
        factor_eta = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "NESTOR Fourier radial branch is invalid")
        if (degree_index < 0 .or. order < 0 .or. scale <= 0.0_dp .or. &
            eta <= 0.0_dp .or. &
            .not. ieee_is_finite(scale) .or. .not. ieee_is_finite(eta)) return
        x = cosh(eta)
        sine_eta = sinh(eta)
        if (use_second_kind) then
            radial_first = toroidal_q_derivative(degree_index, order, x)
            radial_second = toroidal_q_second_derivative(degree_index, order, x)
        else
            radial_first = toroidal_p_derivative(degree_index, order, x)
            radial_second = toroidal_p_second_derivative(degree_index, order, x)
        end if
        if (.not. ieee_is_finite(radial_first) .or. &
            .not. ieee_is_finite(radial_second)) return
        factor = -sine_eta/scale*radial_first
        factor_scale = sine_eta/scale**2*radial_first
        factor_eta = -cosh(eta)/scale*radial_first - &
            sine_eta**2/scale*radial_second
        if (.not. ieee_is_finite(factor) .or. .not. ieee_is_finite(factor_scale) .or. &
            .not. ieee_is_finite(factor_eta)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_factor

    logical function valid_inputs( &
            degree_indices, orders, source_coefficients, scale, eta, output) &
            result(valid)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: source_coefficients(:)
        real(dp), intent(in) :: scale, eta
        complex(dp), intent(in) :: output(:)

        valid = size(degree_indices) > 0 .and. &
            size(orders) == size(degree_indices) .and. &
            size(source_coefficients) == size(degree_indices) .and. &
            size(output) == size(source_coefficients) .and. &
            all(degree_indices >= 0) .and. &
            all(orders >= 0) .and. scale > 0.0_dp .and. eta > 0.0_dp .and. &
            ieee_is_finite(scale) .and. ieee_is_finite(eta) .and. &
            finite_complex(source_coefficients)
    end function valid_inputs

    pure logical function finite_complex(values) result(valid)
        complex(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex

end module fortfem_nestor_fourier_response
