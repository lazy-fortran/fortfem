module fortfem_step_reduction
    !! Smooth actual/predicted reduction diagnostic for nonlinear clients.
    !!
    !! The ratio is intentionally returned without an acceptance decision:
    !! line-search and trust-region policies remain caller-owned and can apply
    !! their own thresholds or topology-event rules.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_step_reduction
    public :: evaluate_step_reduction_jvp
    public :: evaluate_step_reduction_vjp

contains

    subroutine evaluate_step_reduction( &
            base_merit, trial_merit, predicted_reduction, actual_reduction, reduction_ratio, status)
        real(dp), intent(in) :: base_merit, trial_merit, predicted_reduction
        real(dp), intent(out) :: actual_reduction, reduction_ratio
        type(fortsparse_status_t), intent(out) :: status

        actual_reduction = 0.0_dp
        reduction_ratio = 0.0_dp
        if (.not. valid_inputs(base_merit, trial_merit, predicted_reduction)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "step reduction has invalid merit or prediction")
            return
        end if
        actual_reduction = base_merit - trial_merit
        reduction_ratio = actual_reduction/predicted_reduction
        if (.not. ieee_is_finite(actual_reduction) .or. &
            .not. ieee_is_finite(reduction_ratio)) then
            actual_reduction = 0.0_dp
            reduction_ratio = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "step reduction is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_step_reduction

    subroutine evaluate_step_reduction_jvp( &
            base_merit, trial_merit, predicted_reduction, base_merit_dot, trial_merit_dot, &
            predicted_reduction_dot, actual_reduction_dot, reduction_ratio_dot, status)
        real(dp), intent(in) :: base_merit, trial_merit, predicted_reduction
        real(dp), intent(in) :: base_merit_dot, trial_merit_dot, predicted_reduction_dot
        real(dp), intent(out) :: actual_reduction_dot, reduction_ratio_dot
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: actual_reduction

        actual_reduction_dot = 0.0_dp
        reduction_ratio_dot = 0.0_dp
        if (.not. valid_inputs(base_merit, trial_merit, predicted_reduction)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "step reduction JVP has invalid base values")
            return
        end if
        if (.not. valid_directions(base_merit_dot, trial_merit_dot, predicted_reduction_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "step reduction JVP has non-finite directions")
            return
        end if
        actual_reduction = base_merit - trial_merit
        actual_reduction_dot = base_merit_dot - trial_merit_dot
        reduction_ratio_dot = actual_reduction_dot/predicted_reduction - &
            actual_reduction*predicted_reduction_dot/predicted_reduction**2
        if (.not. ieee_is_finite(actual_reduction_dot) .or. &
            .not. ieee_is_finite(reduction_ratio_dot)) then
            actual_reduction_dot = 0.0_dp
            reduction_ratio_dot = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "step reduction JVP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_step_reduction_jvp

    subroutine evaluate_step_reduction_vjp( &
            base_merit, trial_merit, predicted_reduction, actual_reduction_bar, &
            reduction_ratio_bar, base_merit_bar, trial_merit_bar, predicted_reduction_bar, status)
        real(dp), intent(in) :: base_merit, trial_merit, predicted_reduction
        real(dp), intent(in) :: actual_reduction_bar, reduction_ratio_bar
        real(dp), intent(out) :: base_merit_bar, trial_merit_bar, predicted_reduction_bar
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: actual_reduction

        base_merit_bar = 0.0_dp
        trial_merit_bar = 0.0_dp
        predicted_reduction_bar = 0.0_dp
        if (.not. valid_inputs(base_merit, trial_merit, predicted_reduction)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "step reduction VJP has invalid base values")
            return
        end if
        if (.not. valid_directions(actual_reduction_bar, reduction_ratio_bar, 0.0_dp)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "step reduction VJP has non-finite cotangents")
            return
        end if
        actual_reduction = base_merit - trial_merit
        base_merit_bar = actual_reduction_bar + reduction_ratio_bar/predicted_reduction
        trial_merit_bar = -actual_reduction_bar - reduction_ratio_bar/predicted_reduction
        predicted_reduction_bar = -reduction_ratio_bar*actual_reduction/ &
            predicted_reduction**2
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_step_reduction_vjp

    logical function valid_inputs(base_merit, trial_merit, predicted_reduction) result(valid)
        real(dp), intent(in) :: base_merit, trial_merit, predicted_reduction

        valid = ieee_is_finite(base_merit) .and. ieee_is_finite(trial_merit) .and. &
            ieee_is_finite(predicted_reduction) .and. base_merit >= 0.0_dp .and. &
            trial_merit >= 0.0_dp .and. predicted_reduction > 0.0_dp
    end function valid_inputs

    logical function valid_directions(base_dot, trial_dot, predicted_dot) result(valid)
        real(dp), intent(in) :: base_dot, trial_dot, predicted_dot

        valid = ieee_is_finite(base_dot) .and. ieee_is_finite(trial_dot) .and. &
            ieee_is_finite(predicted_dot)
    end function valid_directions

end module fortfem_step_reduction
