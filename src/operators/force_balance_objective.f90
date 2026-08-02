module fortfem_force_balance_objective
    !! Closure-neutral weighted least-squares objective for force residuals.
    !!
    !! The objective is
    !!
    !!   J = 1/2 sum_i w_i r_i^2,
    !!
    !! where a DESC-like or other optimization client supplies the residual
    !! samples and positive quadrature/normalization weights.  No force law,
    !! profile, coordinate convention, or optimizer is selected here.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_force_balance_objective
    public :: evaluate_force_balance_objective_jvp
    public :: evaluate_force_balance_objective_vjp

contains

    subroutine evaluate_force_balance_objective( &
            residual, weights, objective, status)
        real(dp), intent(in) :: residual(:), weights(:)
        real(dp), intent(out) :: objective
        type(fortsparse_status_t), intent(out) :: status

        objective = 0.0_dp
        if (.not. valid_inputs(residual, weights)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "force-balance objective received incompatible samples")
            return
        end if
        objective = 0.5_dp*sum(weights*residual**2)
        if (.not. ieee_is_finite(objective)) then
            objective = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "force-balance objective is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_force_balance_objective

    subroutine evaluate_force_balance_objective_jvp( &
            residual, weights, residual_dot, weights_dot, objective_dot, status)
        real(dp), intent(in) :: residual(:), weights(:), residual_dot(:), weights_dot(:)
        real(dp), intent(out) :: objective_dot
        type(fortsparse_status_t), intent(out) :: status

        objective_dot = 0.0_dp
        if (.not. valid_inputs(residual, weights) .or. &
            size(residual_dot) /= size(residual) .or. size(weights_dot) /= size(weights) .or. &
            .not. all(ieee_is_finite(residual_dot)) .or. &
            .not. all(ieee_is_finite(weights_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "force-balance objective JVP received incompatible samples")
            return
        end if
        objective_dot = sum(weights*residual*residual_dot) + &
            0.5_dp*sum(weights_dot*residual**2)
        if (.not. ieee_is_finite(objective_dot)) then
            objective_dot = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "force-balance objective JVP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_force_balance_objective_jvp

    subroutine evaluate_force_balance_objective_vjp( &
            residual, weights, objective_bar, residual_bar, weights_bar, status)
        real(dp), intent(in) :: residual(:), weights(:), objective_bar
        real(dp), intent(out) :: residual_bar(:), weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        residual_bar = 0.0_dp
        weights_bar = 0.0_dp
        if (.not. valid_inputs(residual, weights) .or. &
            .not. ieee_is_finite(objective_bar) .or. &
            size(residual_bar) /= size(residual) .or. &
            size(weights_bar) /= size(weights)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "force-balance objective VJP received incompatible samples")
            return
        end if
        residual_bar = objective_bar*weights*residual
        weights_bar = 0.5_dp*objective_bar*residual**2
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_force_balance_objective_vjp

    logical function valid_inputs(residual, weights) result(valid)
        real(dp), intent(in) :: residual(:), weights(:)

        valid = size(residual) > 0 .and. size(weights) == size(residual)
        if (.not. valid) return
        valid = all(ieee_is_finite(residual)) .and. &
            all(ieee_is_finite(weights)) .and. all(weights > 0.0_dp)
    end function valid_inputs

end module fortfem_force_balance_objective
