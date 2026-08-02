program test_force_balance_objective
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_force_balance_objective, &
        evaluate_force_balance_objective_jvp, &
        evaluate_force_balance_objective_vjp, &
        evaluate_force_balance_objective_hvp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: residual(3) = [1.0_dp, -2.0_dp, 0.5_dp]
    real(dp), parameter :: weights(3) = [2.0_dp, 0.5_dp, 3.0_dp]
    real(dp), parameter :: residual_dot(3) = [-0.4_dp, 0.7_dp, 0.2_dp]
    real(dp), parameter :: weights_dot(3) = [0.3_dp, -0.1_dp, 0.4_dp]
    real(dp), parameter :: epsilon_fd = 1.0e-6_dp
    real(dp) :: objective, objective_plus, objective_minus, objective_dot
    real(dp) :: residual_bar(3), weights_bar(3), objective_bar
    real(dp) :: residual_bar_dot(3), weights_bar_dot(3)
    real(dp) :: residual_bar_plus(3), weights_bar_plus(3)
    real(dp) :: residual_bar_minus(3), weights_bar_minus(3)
    real(dp) :: lhs, rhs
    logical :: all_passed
    type(fortsparse_status_t) :: status

    all_passed = .true.
    call evaluate_force_balance_objective(residual, weights, objective, status)
    call record_condition(status%code == 0 .and. abs(objective - 2.375_dp) < 1.0e-14_dp, &
        "weighted force objective matches the explicit quadratic oracle")
    call evaluate_force_balance_objective_jvp( &
        residual, weights, residual_dot, weights_dot, objective_dot, status)
    call evaluate_force_balance_objective( &
        residual + epsilon_fd*residual_dot, weights + epsilon_fd*weights_dot, &
        objective_plus, status)
    call evaluate_force_balance_objective( &
        residual - epsilon_fd*residual_dot, weights - epsilon_fd*weights_dot, &
        objective_minus, status)
    call record_condition(status%code == 0 .and. abs(objective_dot - &
        (objective_plus - objective_minus)/(2.0_dp*epsilon_fd)) < 1.0e-8_dp, &
        "weighted force objective JVP matches an independent finite difference")

    objective_bar = -1.7_dp
    call evaluate_force_balance_objective_vjp( &
        residual, weights, objective_bar, residual_bar, weights_bar, status)
    lhs = objective_bar*objective_dot
    rhs = sum(residual_bar*residual_dot) + sum(weights_bar*weights_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-14_dp, &
        "weighted force objective VJP satisfies the real adjoint identity")

    call evaluate_force_balance_objective_hvp( &
        residual, weights, residual_dot, weights_dot, objective_bar, &
        residual_bar_dot, weights_bar_dot, status)
    call evaluate_force_balance_objective_vjp( &
        residual + epsilon_fd*residual_dot, weights + epsilon_fd*weights_dot, &
        objective_bar, residual_bar_plus, weights_bar_plus, status)
    call evaluate_force_balance_objective_vjp( &
        residual - epsilon_fd*residual_dot, weights - epsilon_fd*weights_dot, &
        objective_bar, residual_bar_minus, weights_bar_minus, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(residual_bar_dot - (residual_bar_plus - residual_bar_minus) / &
            (2.0_dp*epsilon_fd))) < 1.0e-8_dp .and. &
        maxval(abs(weights_bar_dot - (weights_bar_plus - weights_bar_minus) / &
            (2.0_dp*epsilon_fd))) < 1.0e-8_dp, &
        "force-balance objective Hessian-vector action matches VJP finite difference")

    call evaluate_force_balance_objective(residual, [-1.0_dp, 1.0_dp, 1.0_dp], &
        objective, status)
    call record_condition(status%code /= 0, &
        "weighted force objective rejects non-positive weights")

    call check_summary("force-balance objective")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_force_balance_objective
