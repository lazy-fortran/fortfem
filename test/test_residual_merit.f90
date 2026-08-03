program test_residual_merit
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_residual_merit, evaluate_residual_merit_jvp, evaluate_residual_merit_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t
    implicit none

    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    logical :: all_passed

    all_passed = .true.
    call check_contract()
    call check_summary("Weighted residual merit")
    if (.not. all_passed) error stop 1

contains

    subroutine check_contract()
        integer, parameter :: residual_size = 4
        real(dp) :: residual(residual_size), weights(residual_size)
        real(dp) :: residual_dot(residual_size), weights_dot(residual_size)
        real(dp) :: residual_plus(residual_size), residual_minus(residual_size)
        real(dp) :: weights_plus(residual_size), weights_minus(residual_size)
        real(dp) :: residual_bar(residual_size), weights_bar(residual_size)
        real(dp) :: merit, merit_dot
        real(dp) :: merit_plus, merit_minus, merit_bar, lhs, rhs
        real(dp) :: expected, expected_dot
        type(fortsparse_status_t) :: status

        residual = [0.4_dp, -0.7_dp, 0.2_dp, 0.9_dp]
        weights = [1.2_dp, 0.8_dp, 1.7_dp, 0.5_dp]
        residual_dot = [-0.3_dp, 0.6_dp, 0.4_dp, -0.2_dp]
        weights_dot = [0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp]
        expected = 0.5_dp*sum(weights*residual**2)
        expected_dot = sum(weights*residual*residual_dot) + &
            0.5_dp*sum(weights_dot*residual**2)

        call evaluate_residual_merit(residual, weights, merit, status)
        call record(status%code == FORTSPARSE_OK .and. abs(merit - expected) < 1.0e-14_dp, &
            "residual merit value oracle")

        call evaluate_residual_merit_jvp( &
            residual, weights, residual_dot, weights_dot, merit_dot, status)
        call evaluate_residual_merit( &
            residual + epsilon_fd*residual_dot, weights + epsilon_fd*weights_dot, &
            merit_plus, status)
        call evaluate_residual_merit( &
            residual - epsilon_fd*residual_dot, weights - epsilon_fd*weights_dot, &
            merit_minus, status)
        call record(status%code == FORTSPARSE_OK .and. abs(merit_dot - expected_dot) < 1.0e-14_dp .and. &
            abs(merit_dot - (merit_plus - merit_minus)/(2.0_dp*epsilon_fd)) < 2.0e-8_dp, &
            "residual merit JVP oracle")

        residual_bar = [0.2_dp, -0.8_dp, 0.5_dp, 0.6_dp]
        merit_bar = -0.7_dp
        call evaluate_residual_merit_vjp( &
            residual, weights, merit_bar, residual_bar, weights_bar, status)
        lhs = merit_bar*merit_dot
        rhs = dot_product(residual_bar, residual_dot) + dot_product(weights_bar, weights_dot)
        call record(status%code == FORTSPARSE_OK .and. abs(lhs - rhs) < 1.0e-13_dp, &
            "residual merit VJP dot-product oracle")

        call evaluate_residual_merit(residual, [1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp], merit, status)
        call record(status%code == FORTSPARSE_INVALID_MATRIX, &
            "non-positive merit weight is rejected")
    end subroutine check_contract

    subroutine record(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        call check_condition(condition, message)
        all_passed = all_passed .and. condition
    end subroutine record

end program test_residual_merit
