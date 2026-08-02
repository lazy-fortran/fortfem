module fortfem_toroidal_diagnostic_hooks
    !! Geometry-neutral diagnostics on caller-supplied toroidal samples.
    !!
    !! The routines in this module intentionally do not construct flux
    !! coordinates, a Boozer map, an equilibrium, or a near-axis expansion.
    !! They reduce data which an application has already sampled on a fixed
    !! topology.  This keeps the contracts usable by nested-surface, Fourier,
    !! IGA, and external equilibrium clients without importing their readers
    !! or constitutive physics.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    type, public :: near_axis_diagnostic_metadata_t
        !! Fixed-topology summary of supplied near-axis coefficients.
        !!
        !! A term with radial power p and poloidal mode m is axis-regular for
        !! this diagnostic when p >= abs(m) and p-abs(m) is even.  This is a
        !! metadata check only; it does not modify or reject an otherwise
        !! usable supplied expansion.  The axis value is the sum of all
        !! coefficients with p=0,m=0.
        character(len=40) :: schema_version = "fortfem-near-axis-diagnostic-1"
        integer :: coefficient_count = 0
        integer :: regular_count = 0
        integer :: irregular_count = 0
        integer :: axis_mode_count = 0
        integer :: max_radial_power = 0
        real(dp) :: coefficient_l2_norm = 0.0_dp
        real(dp) :: axis_value_norm = 0.0_dp
        complex(dp) :: axis_value = cmplx(0.0_dp, 0.0_dp, dp)
        logical :: axis_regular = .false.
    end type near_axis_diagnostic_metadata_t

    public :: evaluate_boozer_like_rotational_transform
    public :: evaluate_boozer_like_rotational_transform_jvp
    public :: evaluate_boozer_like_rotational_transform_vjp
    public :: evaluate_near_axis_diagnostic_metadata
    public :: evaluate_near_axis_diagnostic_metadata_jvp
    public :: evaluate_near_axis_diagnostic_metadata_vjp

contains

    subroutine evaluate_boozer_like_rotational_transform( &
            poloidal_rate, toroidal_rate, weights, rotational_transform, &
            poloidal_integral, toroidal_integral, status)
        !! Form a neutral weighted rotational-transform reduction.
        !!
        !! The caller supplies rates that have whatever field/coordinate
        !! convention it has selected.  The result is
        !!
        !!   iota = sum(w * poloidal_rate) / sum(w * toroidal_rate).
        !!
        !! It is therefore Boozer-like rather than a claim to construct a
        !! Boozer coordinate transformation.  The two integrals are returned
        !! to make normalization and orientation choices inspectable.
        real(dp), intent(in) :: poloidal_rate(:), toroidal_rate(:), weights(:)
        real(dp), intent(out) :: rotational_transform
        real(dp), intent(out) :: poloidal_integral, toroidal_integral
        type(fortsparse_status_t), intent(out) :: status

        rotational_transform = 0.0_dp
        poloidal_integral = 0.0_dp
        toroidal_integral = 0.0_dp
        if (.not. valid_rate_inputs(poloidal_rate, toroidal_rate, weights)) then
            call invalid_status(status, &
                "rotational-transform rates and weights are incompatible")
            return
        end if
        poloidal_integral = sum(weights*poloidal_rate)
        toroidal_integral = sum(weights*toroidal_rate)
        if (.not. ieee_is_finite(poloidal_integral) .or. &
            .not. ieee_is_finite(toroidal_integral) .or. &
            abs(toroidal_integral) <= tiny(1.0_dp)) then
            poloidal_integral = 0.0_dp
            toroidal_integral = 0.0_dp
            call invalid_status(status, &
                "rotational-transform toroidal integral is singular")
            return
        end if
        rotational_transform = poloidal_integral/toroidal_integral
        if (.not. ieee_is_finite(rotational_transform)) then
            rotational_transform = 0.0_dp
            poloidal_integral = 0.0_dp
            toroidal_integral = 0.0_dp
            call invalid_status(status, "rotational-transform reduction is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_boozer_like_rotational_transform

    subroutine evaluate_boozer_like_rotational_transform_jvp( &
            poloidal_rate, toroidal_rate, weights, poloidal_rate_dot, &
            toroidal_rate_dot, weights_dot, rotational_transform_dot, &
            poloidal_integral_dot, toroidal_integral_dot, status)
        !! Fixed-topology JVP of the supplied-rate rotational reduction.
        real(dp), intent(in) :: poloidal_rate(:), toroidal_rate(:), weights(:)
        real(dp), intent(in) :: poloidal_rate_dot(:), toroidal_rate_dot(:), weights_dot(:)
        real(dp), intent(out) :: rotational_transform_dot
        real(dp), intent(out) :: poloidal_integral_dot, toroidal_integral_dot
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: rotational_transform, poloidal_integral, toroidal_integral

        rotational_transform_dot = 0.0_dp
        poloidal_integral_dot = 0.0_dp
        toroidal_integral_dot = 0.0_dp
        call evaluate_boozer_like_rotational_transform( &
            poloidal_rate, toroidal_rate, weights, rotational_transform, &
            poloidal_integral, toroidal_integral, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_rate_directions( &
                poloidal_rate, toroidal_rate, weights, poloidal_rate_dot, &
                toroidal_rate_dot, weights_dot)) then
            call invalid_status(status, &
                "rotational-transform JVP directions are incompatible")
            return
        end if
        poloidal_integral_dot = sum(weights_dot*poloidal_rate + &
            weights*poloidal_rate_dot)
        toroidal_integral_dot = sum(weights_dot*toroidal_rate + &
            weights*toroidal_rate_dot)
        rotational_transform_dot = (poloidal_integral_dot - rotational_transform* &
            toroidal_integral_dot)/toroidal_integral
        if (.not. ieee_is_finite(rotational_transform_dot) .or. &
            .not. ieee_is_finite(poloidal_integral_dot) .or. &
            .not. ieee_is_finite(toroidal_integral_dot)) then
            rotational_transform_dot = 0.0_dp
            poloidal_integral_dot = 0.0_dp
            toroidal_integral_dot = 0.0_dp
            call invalid_status(status, "rotational-transform JVP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_boozer_like_rotational_transform_jvp

    subroutine evaluate_boozer_like_rotational_transform_vjp( &
            poloidal_rate, toroidal_rate, weights, rotational_transform_bar, &
            poloidal_rate_bar, toroidal_rate_bar, weights_bar, status)
        !! Real adjoint of the scalar rotational-transform reduction.
        real(dp), intent(in) :: poloidal_rate(:), toroidal_rate(:), weights(:)
        real(dp), intent(in) :: rotational_transform_bar
        real(dp), intent(out) :: poloidal_rate_bar(:), toroidal_rate_bar(:), weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: rotational_transform, poloidal_integral, toroidal_integral
        real(dp) :: numerator_bar, denominator_bar

        poloidal_rate_bar = 0.0_dp
        toroidal_rate_bar = 0.0_dp
        weights_bar = 0.0_dp
        call evaluate_boozer_like_rotational_transform( &
            poloidal_rate, toroidal_rate, weights, rotational_transform, &
            poloidal_integral, toroidal_integral, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. ieee_is_finite(rotational_transform_bar) .or. &
            size(poloidal_rate_bar) /= size(poloidal_rate) .or. &
            size(toroidal_rate_bar) /= size(toroidal_rate) .or. &
            size(weights_bar) /= size(weights)) then
            call invalid_status(status, &
                "rotational-transform VJP outputs or cotangent are invalid")
            return
        end if
        numerator_bar = rotational_transform_bar/toroidal_integral
        denominator_bar = -rotational_transform_bar*poloidal_integral/ &
            toroidal_integral**2
        poloidal_rate_bar = numerator_bar*weights
        toroidal_rate_bar = denominator_bar*weights
        weights_bar = numerator_bar*poloidal_rate + denominator_bar*toroidal_rate
        if (.not. all(ieee_is_finite(poloidal_rate_bar)) .or. &
            .not. all(ieee_is_finite(toroidal_rate_bar)) .or. &
            .not. all(ieee_is_finite(weights_bar))) then
            poloidal_rate_bar = 0.0_dp
            toroidal_rate_bar = 0.0_dp
            weights_bar = 0.0_dp
            call invalid_status(status, "rotational-transform VJP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_boozer_like_rotational_transform_vjp

    subroutine evaluate_near_axis_diagnostic_metadata( &
            radial_powers, poloidal_modes, coefficients, metadata, status)
        !! Summarize caller-supplied near-axis/Fourier-Zernike coefficients.
        integer, intent(in) :: radial_powers(:), poloidal_modes(:)
        complex(dp), intent(in) :: coefficients(:)
        type(near_axis_diagnostic_metadata_t), intent(out) :: metadata
        type(fortsparse_status_t), intent(out) :: status
        integer :: term

        metadata = near_axis_diagnostic_metadata_t()
        if (.not. valid_coefficient_inputs( &
                radial_powers, poloidal_modes, coefficients)) then
            call invalid_status(status, &
                "near-axis coefficients and mode metadata are incompatible")
            return
        end if
        metadata%coefficient_count = size(coefficients)
        metadata%max_radial_power = maxval(radial_powers)
        do term = 1, size(coefficients)
            metadata%coefficient_l2_norm = metadata%coefficient_l2_norm + &
                abs(coefficients(term))**2
            if (radial_powers(term) == 0 .and. poloidal_modes(term) == 0) then
                metadata%axis_mode_count = metadata%axis_mode_count + 1
                metadata%axis_value = metadata%axis_value + coefficients(term)
            end if
            if (is_axis_regular_term( &
                    radial_powers(term), poloidal_modes(term))) then
                metadata%regular_count = metadata%regular_count + 1
            else
                metadata%irregular_count = metadata%irregular_count + 1
            end if
        end do
        metadata%coefficient_l2_norm = sqrt(metadata%coefficient_l2_norm)
        metadata%axis_value_norm = abs(metadata%axis_value)
        metadata%axis_regular = metadata%irregular_count == 0
        if (.not. ieee_is_finite(metadata%coefficient_l2_norm) .or. &
            .not. ieee_is_finite(metadata%axis_value_norm)) then
            metadata = near_axis_diagnostic_metadata_t()
            call invalid_status(status, "near-axis coefficient metadata is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_near_axis_diagnostic_metadata

    subroutine evaluate_near_axis_diagnostic_metadata_jvp( &
            radial_powers, poloidal_modes, coefficients, coefficients_dot, &
            coefficient_l2_norm_dot, axis_value_dot, axis_value_norm_dot, status)
        !! Derivatives of the smooth scalar metadata fields.
        integer, intent(in) :: radial_powers(:), poloidal_modes(:)
        complex(dp), intent(in) :: coefficients(:), coefficients_dot(:)
        real(dp), intent(out) :: coefficient_l2_norm_dot, axis_value_norm_dot
        complex(dp), intent(out) :: axis_value_dot
        type(fortsparse_status_t), intent(out) :: status
        type(near_axis_diagnostic_metadata_t) :: metadata

        coefficient_l2_norm_dot = 0.0_dp
        axis_value_dot = cmplx(0.0_dp, 0.0_dp, dp)
        axis_value_norm_dot = 0.0_dp
        call evaluate_near_axis_diagnostic_metadata( &
            radial_powers, poloidal_modes, coefficients, metadata, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(coefficients_dot) /= size(coefficients) .or. &
            .not. finite_complex(coefficients_dot)) then
            call invalid_status(status, &
                "near-axis coefficient JVP direction is incompatible")
            return
        end if
        coefficient_l2_norm_dot = real(sum(conjg(coefficients)*coefficients_dot), dp) / &
            max(metadata%coefficient_l2_norm, tiny(1.0_dp))
        do concurrent (integer :: term = 1:size(coefficients))
            if (radial_powers(term) == 0 .and. poloidal_modes(term) == 0) then
                axis_value_dot = axis_value_dot + coefficients_dot(term)
            end if
        end do
        if (metadata%axis_value_norm > tiny(1.0_dp)) then
            axis_value_norm_dot = real(conjg(metadata%axis_value)*axis_value_dot, dp) / &
                metadata%axis_value_norm
        end if
        if (.not. ieee_is_finite(coefficient_l2_norm_dot) .or. &
            .not. ieee_is_finite(axis_value_norm_dot) .or. &
            .not. finite_complex([axis_value_dot])) then
            coefficient_l2_norm_dot = 0.0_dp
            axis_value_dot = cmplx(0.0_dp, 0.0_dp, dp)
            axis_value_norm_dot = 0.0_dp
            call invalid_status(status, "near-axis coefficient JVP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_near_axis_diagnostic_metadata_jvp

    subroutine evaluate_near_axis_diagnostic_metadata_vjp( &
            radial_powers, poloidal_modes, coefficients, coefficient_l2_norm_bar, &
            axis_value_bar, axis_value_norm_bar, coefficients_bar, status)
        !! Real-part complex adjoint for the smooth metadata fields.
        integer, intent(in) :: radial_powers(:), poloidal_modes(:)
        complex(dp), intent(in) :: coefficients(:), axis_value_bar
        real(dp), intent(in) :: coefficient_l2_norm_bar, axis_value_norm_bar
        complex(dp), intent(out) :: coefficients_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        type(near_axis_diagnostic_metadata_t) :: metadata
        integer :: term
        complex(dp) :: axis_norm_gradient

        coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        call evaluate_near_axis_diagnostic_metadata( &
            radial_powers, poloidal_modes, coefficients, metadata, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(coefficients_bar) /= size(coefficients) .or. &
            .not. ieee_is_finite(coefficient_l2_norm_bar) .or. &
            .not. ieee_is_finite(axis_value_norm_bar) .or. &
            .not. finite_complex([axis_value_bar])) then
            call invalid_status(status, "near-axis coefficient VJP inputs are invalid")
            return
        end if
        if (metadata%coefficient_l2_norm > tiny(1.0_dp)) then
            coefficients_bar = coefficient_l2_norm_bar*coefficients / &
                metadata%coefficient_l2_norm
        end if
        if (metadata%axis_value_norm > tiny(1.0_dp)) then
            axis_norm_gradient = axis_value_norm_bar*metadata%axis_value / &
                metadata%axis_value_norm
        else
            axis_norm_gradient = cmplx(0.0_dp, 0.0_dp, dp)
        end if
        do term = 1, size(coefficients)
            if (radial_powers(term) == 0 .and. poloidal_modes(term) == 0) then
                coefficients_bar(term) = coefficients_bar(term) + axis_value_bar + &
                    axis_norm_gradient
            end if
        end do
        if (.not. finite_complex(coefficients_bar)) then
            coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
            call invalid_status(status, "near-axis coefficient VJP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_near_axis_diagnostic_metadata_vjp

    logical function valid_rate_inputs( &
            poloidal_rate, toroidal_rate, weights) result(valid)
        real(dp), intent(in) :: poloidal_rate(:), toroidal_rate(:), weights(:)

        valid = size(poloidal_rate) > 0 .and. &
            size(toroidal_rate) == size(poloidal_rate) .and. &
            size(weights) == size(poloidal_rate) .and. &
            all(ieee_is_finite(poloidal_rate)) .and. &
            all(ieee_is_finite(toroidal_rate)) .and. &
            all(ieee_is_finite(weights)) .and. all(weights > 0.0_dp)
    end function valid_rate_inputs

    logical function valid_rate_directions( &
            poloidal_rate, toroidal_rate, weights, poloidal_rate_dot, &
            toroidal_rate_dot, weights_dot) result(valid)
        real(dp), intent(in) :: poloidal_rate(:), toroidal_rate(:), weights(:)
        real(dp), intent(in) :: poloidal_rate_dot(:), toroidal_rate_dot(:), weights_dot(:)

        valid = valid_rate_inputs(poloidal_rate, toroidal_rate, weights) .and. &
            size(poloidal_rate_dot) == size(poloidal_rate) .and. &
            size(toroidal_rate_dot) == size(toroidal_rate) .and. &
            size(weights_dot) == size(weights) .and. &
            all(ieee_is_finite(poloidal_rate_dot)) .and. &
            all(ieee_is_finite(toroidal_rate_dot)) .and. &
            all(ieee_is_finite(weights_dot))
    end function valid_rate_directions

    logical function valid_coefficient_inputs( &
            radial_powers, poloidal_modes, coefficients) result(valid)
        integer, intent(in) :: radial_powers(:), poloidal_modes(:)
        complex(dp), intent(in) :: coefficients(:)

        valid = size(coefficients) > 0 .and. &
            size(radial_powers) == size(coefficients) .and. &
            size(poloidal_modes) == size(coefficients) .and. &
            all(radial_powers >= 0) .and. finite_complex(coefficients)
    end function valid_coefficient_inputs

    pure logical function is_axis_regular_term(radial_power, poloidal_mode) result(valid)
        integer, intent(in) :: radial_power, poloidal_mode

        valid = radial_power >= abs(poloidal_mode) .and. &
            mod(radial_power - abs(poloidal_mode), 2) == 0
    end function is_axis_regular_term

    pure logical function finite_complex(values) result(valid)
        complex(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex

    subroutine invalid_status(status, message)
        type(fortsparse_status_t), intent(out) :: status
        character(len=*), intent(in) :: message

        call status_set(status, FORTSPARSE_INVALID_MATRIX, message)
    end subroutine invalid_status

end module fortfem_toroidal_diagnostic_hooks
