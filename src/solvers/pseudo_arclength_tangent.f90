module fortfem_pseudo_arclength_tangent
    !! Fixed-topology normalization of a state/parameter continuation tangent.
    !!
    !! The tangent is treated as one Euclidean vector split into state and
    !! parameter blocks.  This primitive supplies the scale used by a
    !! pseudo-arclength residual without choosing a continuation policy.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: normalize_pseudo_arclength_tangent
    public :: normalize_pseudo_arclength_tangent_jvp
    public :: normalize_pseudo_arclength_tangent_vjp

contains

    subroutine normalize_pseudo_arclength_tangent( &
            tangent_state, tangent_parameter, normalized_state, normalized_parameter, &
            norm, status)
        real(dp), intent(in) :: tangent_state(:), tangent_parameter(:)
        real(dp), intent(out) :: normalized_state(:), normalized_parameter(:), norm
        type(fortsparse_status_t), intent(out) :: status

        normalized_state = 0.0_dp
        normalized_parameter = 0.0_dp
        norm = 0.0_dp
        if (.not. valid_base_inputs( &
            tangent_state, tangent_parameter, normalized_state, normalized_parameter)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-arclength tangent has incompatible inputs")
            return
        end if
        norm = sqrt(dot_product(tangent_state, tangent_state) + &
            dot_product(tangent_parameter, tangent_parameter))
        if (.not. ieee_is_finite(norm) .or. norm <= tiny(1.0_dp)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-arclength tangent has zero or non-finite norm")
            return
        end if
        normalized_state = tangent_state/norm
        normalized_parameter = tangent_parameter/norm
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine normalize_pseudo_arclength_tangent

    subroutine normalize_pseudo_arclength_tangent_jvp( &
            tangent_state, tangent_parameter, tangent_state_dot, tangent_parameter_dot, &
            normalized_state_dot, normalized_parameter_dot, norm_dot, status)
        real(dp), intent(in) :: tangent_state(:), tangent_parameter(:)
        real(dp), intent(in) :: tangent_state_dot(:), tangent_parameter_dot(:)
        real(dp), intent(out) :: normalized_state_dot(:), normalized_parameter_dot(:), norm_dot
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: norm, tangent_projection

        normalized_state_dot = 0.0_dp
        normalized_parameter_dot = 0.0_dp
        norm_dot = 0.0_dp
        if (.not. valid_base_inputs( &
            tangent_state, tangent_parameter, normalized_state_dot, normalized_parameter_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-arclength tangent JVP has incompatible base inputs")
            return
        end if
        if (.not. valid_direction_inputs( &
            tangent_state, tangent_parameter, tangent_state_dot, tangent_parameter_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-arclength tangent JVP has incompatible directions")
            return
        end if
        norm = sqrt(dot_product(tangent_state, tangent_state) + &
            dot_product(tangent_parameter, tangent_parameter))
        if (.not. ieee_is_finite(norm) .or. norm <= tiny(1.0_dp)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-arclength tangent JVP has zero or non-finite norm")
            return
        end if
        tangent_projection = dot_product(tangent_state, tangent_state_dot) + &
            dot_product(tangent_parameter, tangent_parameter_dot)
        norm_dot = tangent_projection/norm
        normalized_state_dot = (tangent_state_dot - tangent_state*norm_dot/norm)/norm
        normalized_parameter_dot = (tangent_parameter_dot - tangent_parameter*norm_dot/norm)/norm
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine normalize_pseudo_arclength_tangent_jvp

    subroutine normalize_pseudo_arclength_tangent_vjp( &
            tangent_state, tangent_parameter, normalized_state_bar, normalized_parameter_bar, &
            norm_bar, tangent_state_bar, tangent_parameter_bar, status)
        real(dp), intent(in) :: tangent_state(:), tangent_parameter(:)
        real(dp), intent(in) :: normalized_state_bar(:), normalized_parameter_bar(:), norm_bar
        real(dp), intent(out) :: tangent_state_bar(:), tangent_parameter_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: norm, cotangent_projection

        tangent_state_bar = 0.0_dp
        tangent_parameter_bar = 0.0_dp
        if (.not. valid_base_inputs( &
            tangent_state, tangent_parameter, tangent_state_bar, tangent_parameter_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-arclength tangent VJP has incompatible base inputs")
            return
        end if
        if (.not. valid_cotangent_inputs( &
            tangent_state, tangent_parameter, normalized_state_bar, normalized_parameter_bar, norm_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-arclength tangent VJP has incompatible cotangents")
            return
        end if
        norm = sqrt(dot_product(tangent_state, tangent_state) + &
            dot_product(tangent_parameter, tangent_parameter))
        if (.not. ieee_is_finite(norm) .or. norm <= tiny(1.0_dp)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "pseudo-arclength tangent VJP has zero or non-finite norm")
            return
        end if
        cotangent_projection = (dot_product(normalized_state_bar, tangent_state) + &
            dot_product(normalized_parameter_bar, tangent_parameter))/norm
        tangent_state_bar = normalized_state_bar/norm - &
            tangent_state*cotangent_projection/(norm*norm) + norm_bar*tangent_state/norm
        tangent_parameter_bar = normalized_parameter_bar/norm - &
            tangent_parameter*cotangent_projection/(norm*norm) + norm_bar*tangent_parameter/norm
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine normalize_pseudo_arclength_tangent_vjp

    logical function valid_base_inputs( &
            tangent_state, tangent_parameter, target_state, target_parameter) result(valid)
        real(dp), intent(in) :: tangent_state(:), tangent_parameter(:)
        real(dp), intent(in) :: target_state(:), target_parameter(:)

        valid = size(tangent_state) > 0 .and. size(tangent_parameter) > 0
        if (.not. valid) return
        if (size(target_state) /= size(tangent_state) .or. &
            size(target_parameter) /= size(tangent_parameter)) then
            valid = .false.
            return
        end if
        valid = all(ieee_is_finite(tangent_state)) .and. &
            all(ieee_is_finite(tangent_parameter)) .and. &
            all(ieee_is_finite(target_state)) .and. all(ieee_is_finite(target_parameter))
    end function valid_base_inputs

    logical function valid_direction_inputs( &
            tangent_state, tangent_parameter, tangent_state_dot, tangent_parameter_dot) result(valid)
        real(dp), intent(in) :: tangent_state(:), tangent_parameter(:)
        real(dp), intent(in) :: tangent_state_dot(:), tangent_parameter_dot(:)

        valid = size(tangent_state_dot) == size(tangent_state) .and. &
            size(tangent_parameter_dot) == size(tangent_parameter)
        if (.not. valid) return
        valid = all(ieee_is_finite(tangent_state_dot)) .and. &
            all(ieee_is_finite(tangent_parameter_dot))
    end function valid_direction_inputs

    logical function valid_cotangent_inputs( &
            tangent_state, tangent_parameter, normalized_state_bar, &
            normalized_parameter_bar, norm_bar) result(valid)
        real(dp), intent(in) :: tangent_state(:), tangent_parameter(:)
        real(dp), intent(in) :: normalized_state_bar(:), normalized_parameter_bar(:), norm_bar

        valid = size(normalized_state_bar) == size(tangent_state) .and. &
            size(normalized_parameter_bar) == size(tangent_parameter)
        if (.not. valid) return
        valid = ieee_is_finite(norm_bar) .and. all(ieee_is_finite(normalized_state_bar)) .and. &
            all(ieee_is_finite(normalized_parameter_bar))
    end function valid_cotangent_inputs

end module fortfem_pseudo_arclength_tangent
