program test_step_reduction
    use check, only: check_condition, check_summary
    use fortfem_step_reduction, only: &
        evaluate_step_reduction, evaluate_step_reduction_jvp, evaluate_step_reduction_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t
    implicit none

    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    logical :: all_passed

    all_passed = .true.
    call check_contract()
    call check_summary("Step reduction diagnostic")
    if (.not. all_passed) error stop 1

contains

    subroutine check_contract()
        real(dp) :: base_merit, trial_merit, predicted_reduction
        real(dp) :: base_merit_dot, trial_merit_dot, predicted_reduction_dot
        real(dp) :: actual_reduction, reduction_ratio, actual_reduction_dot, reduction_ratio_dot
        real(dp) :: actual_reduction_plus, actual_reduction_minus
        real(dp) :: reduction_ratio_plus, reduction_ratio_minus
        real(dp) :: actual_reduction_bar, reduction_ratio_bar
        real(dp) :: base_merit_bar, trial_merit_bar, predicted_reduction_bar
        real(dp) :: expected_actual, expected_ratio, expected_actual_dot, expected_ratio_dot
        real(dp) :: lhs, rhs
        type(fortsparse_status_t) :: status

        base_merit = 2.4_dp
        trial_merit = 1.1_dp
        predicted_reduction = 1.6_dp
        base_merit_dot = -0.3_dp
        trial_merit_dot = 0.5_dp
        predicted_reduction_dot = 0.2_dp
        expected_actual = base_merit - trial_merit
        expected_ratio = expected_actual/predicted_reduction
        expected_actual_dot = base_merit_dot - trial_merit_dot
        expected_ratio_dot = expected_actual_dot/predicted_reduction - &
            expected_actual*predicted_reduction_dot/predicted_reduction**2

        call evaluate_step_reduction( &
            base_merit, trial_merit, predicted_reduction, actual_reduction, reduction_ratio, status)
        call record(status%code == FORTSPARSE_OK .and. &
            abs(actual_reduction - expected_actual) < 1.0e-14_dp .and. &
            abs(reduction_ratio - expected_ratio) < 1.0e-14_dp, &
            "step reduction value oracle")

        call evaluate_step_reduction_jvp( &
            base_merit, trial_merit, predicted_reduction, base_merit_dot, trial_merit_dot, &
            predicted_reduction_dot, actual_reduction_dot, reduction_ratio_dot, status)
        call evaluate_step_reduction( &
            base_merit + epsilon_fd*base_merit_dot, trial_merit + epsilon_fd*trial_merit_dot, &
            predicted_reduction + epsilon_fd*predicted_reduction_dot, actual_reduction_plus, &
            reduction_ratio_plus, status)
        call evaluate_step_reduction( &
            base_merit - epsilon_fd*base_merit_dot, trial_merit - epsilon_fd*trial_merit_dot, &
            predicted_reduction - epsilon_fd*predicted_reduction_dot, actual_reduction_minus, &
            reduction_ratio_minus, status)
        call record(status%code == FORTSPARSE_OK .and. &
            abs(actual_reduction_dot - expected_actual_dot) < 1.0e-14_dp .and. &
            abs(reduction_ratio_dot - expected_ratio_dot) < 1.0e-14_dp .and. &
            abs(actual_reduction_dot - (actual_reduction_plus - actual_reduction_minus)/ &
            (2.0_dp*epsilon_fd)) < 2.0e-8_dp .and. &
            abs(reduction_ratio_dot - (reduction_ratio_plus - reduction_ratio_minus)/ &
            (2.0_dp*epsilon_fd)) < 2.0e-8_dp, "step reduction JVP oracle")

        actual_reduction_bar = 0.4_dp
        reduction_ratio_bar = -0.7_dp
        call evaluate_step_reduction_vjp( &
            base_merit, trial_merit, predicted_reduction, actual_reduction_bar, reduction_ratio_bar, &
            base_merit_bar, trial_merit_bar, predicted_reduction_bar, status)
        lhs = actual_reduction_bar*actual_reduction_dot + reduction_ratio_bar*reduction_ratio_dot
        rhs = base_merit_bar*base_merit_dot + trial_merit_bar*trial_merit_dot + &
            predicted_reduction_bar*predicted_reduction_dot
        call record(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < 1.0e-13_dp, &
            "step reduction VJP dot-product oracle")

        call evaluate_step_reduction( &
            base_merit, trial_merit, 0.0_dp, actual_reduction, reduction_ratio, status)
        call record(status%code == FORTSPARSE_INVALID_MATRIX, &
            "non-positive predicted reduction is rejected")
    end subroutine check_contract

    subroutine record(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        call check_condition(condition, message)
        all_passed = all_passed .and. condition
    end subroutine record

end program test_step_reduction
