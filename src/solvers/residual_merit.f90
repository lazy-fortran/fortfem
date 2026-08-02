module fortfem_residual_merit
    !! Weighted least-squares merit for line-search and trust-region clients.
    !!
    !! The merit is a scalar diagnostic only; this module does not accept or
    !! reject a step and does not own a nonlinear solver or stopping policy.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_residual_merit
    public :: evaluate_residual_merit_jvp
    public :: evaluate_residual_merit_vjp

contains

    subroutine evaluate_residual_merit(residual, weights, merit, status)
        real(dp), intent(in) :: residual(:), weights(:)
        real(dp), intent(out) :: merit
        type(fortsparse_status_t), intent(out) :: status

        merit = 0.0_dp
        if (.not. valid_inputs(residual, weights)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "residual merit has incompatible or non-positive inputs")
            return
        end if
        merit = 0.5_dp*sum(weights*residual**2)
        if (.not. ieee_is_finite(merit)) then
            merit = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "residual merit is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_residual_merit

    subroutine evaluate_residual_merit_jvp( &
            residual, weights, residual_dot, weights_dot, merit_dot, status)
        real(dp), intent(in) :: residual(:), weights(:)
        real(dp), intent(in) :: residual_dot(:), weights_dot(:)
        real(dp), intent(out) :: merit_dot
        type(fortsparse_status_t), intent(out) :: status

        merit_dot = 0.0_dp
        if (.not. valid_inputs(residual, weights)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "residual merit JVP has incompatible base inputs")
            return
        end if
        if (.not. valid_directions(residual, weights, residual_dot, weights_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "residual merit JVP has incompatible directions")
            return
        end if
        merit_dot = sum(weights*residual*residual_dot) + &
            0.5_dp*sum(weights_dot*residual**2)
        if (.not. ieee_is_finite(merit_dot)) then
            merit_dot = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "residual merit JVP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_residual_merit_jvp

    subroutine evaluate_residual_merit_vjp( &
            residual, weights, merit_bar, residual_bar, weights_bar, status)
        real(dp), intent(in) :: residual(:), weights(:), merit_bar
        real(dp), intent(out) :: residual_bar(:), weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        residual_bar = 0.0_dp
        weights_bar = 0.0_dp
        if (.not. valid_inputs(residual, weights)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "residual merit VJP has incompatible base inputs")
            return
        end if
        if (.not. ieee_is_finite(merit_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "residual merit VJP has a non-finite cotangent")
            return
        end if
        if (size(residual_bar) /= size(residual) .or. &
            size(weights_bar) /= size(weights)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "residual merit VJP has incompatible cotangent shapes")
            return
        end if
        residual_bar = merit_bar*weights*residual
        weights_bar = 0.5_dp*merit_bar*residual**2
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_residual_merit_vjp

    logical function valid_inputs(residual, weights) result(valid)
        real(dp), intent(in) :: residual(:), weights(:)

        valid = size(residual) > 0 .and. size(weights) == size(residual)
        if (.not. valid) return
        if (.not. all(ieee_is_finite(residual)) .or. &
            .not. all(ieee_is_finite(weights)) .or. any(weights <= 0.0_dp)) then
            valid = .false.
        end if
    end function valid_inputs

    logical function valid_directions( &
            residual, weights, residual_dot, weights_dot) result(valid)
        real(dp), intent(in) :: residual(:), weights(:)
        real(dp), intent(in) :: residual_dot(:), weights_dot(:)

        valid = size(residual_dot) == size(residual) .and. &
            size(weights_dot) == size(weights)
        if (.not. valid) return
        valid = all(ieee_is_finite(residual_dot)) .and. &
            all(ieee_is_finite(weights_dot))
    end function valid_directions

end module fortfem_residual_merit
