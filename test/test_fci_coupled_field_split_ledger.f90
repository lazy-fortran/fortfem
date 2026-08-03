program test_fci_coupled_field_split_ledger
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_fci_coupled_field_split_ledger, &
        evaluate_fci_coupled_field_split_ledger_jvp, &
        evaluate_fci_coupled_field_split_ledger_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: step = 1.0e-5_dp
    real(dp), parameter :: residual(3) = [1.0_dp, -2.0_dp, 0.5_dp]
    real(dp), parameter :: parallel_action(3) = [0.5_dp, -0.5_dp, 0.2_dp]
    real(dp), parameter :: plane_action(3) = [0.2_dp, -0.8_dp, 0.6_dp]
    real(dp), parameter :: retained_actions(3, 2) = reshape([ &
        1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, 0.5_dp], [3, 2])
    real(dp), parameter :: weights(4) = [0.4_dp, 0.7_dp, 1.1_dp, 0.3_dp]
    real(dp), parameter :: residual_dot(3) = [0.1_dp, 0.2_dp, -0.3_dp]
    real(dp), parameter :: parallel_action_dot(3) = [-0.2_dp, 0.1_dp, 0.05_dp]
    real(dp), parameter :: plane_action_dot(3) = [0.3_dp, -0.1_dp, 0.2_dp]
    real(dp), parameter :: retained_actions_dot(3, 2) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp, 0.2_dp, 0.1_dp], [3, 2])
    real(dp), parameter :: weights_dot(4) = [0.05_dp, -0.03_dp, 0.02_dp, 0.04_dp]
    real(dp), parameter :: correction_bar(3) = [0.7_dp, -0.4_dp, 0.2_dp]
    real(dp), parameter :: ledger_bar(4) = [0.1_dp, -0.5_dp, 0.6_dp, 0.3_dp]
    real(dp), parameter :: total_bar = -0.25_dp
    real(dp) :: correction(3), ledger(4), total
    real(dp) :: correction_dot(3), ledger_dot(4), total_dot
    real(dp) :: correction_plus(3), correction_minus(3)
    real(dp) :: ledger_plus(4), ledger_minus(4), total_plus, total_minus
    real(dp) :: residual_bar(3), parallel_action_bar(3), plane_action_bar(3)
    real(dp) :: retained_actions_bar(3, 2), weights_bar(4)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status

    call evaluate_fci_coupled_field_split_ledger( &
        residual, parallel_action, plane_action, retained_actions, weights, &
        correction, ledger, total, status)
    call check_condition(status%code == 0, &
        "FCI coupled split accepts positive caller-owned actions")
    call check_condition(maxval(abs(correction - [1.44_dp, -1.06_dp, 1.75_dp])) &
        < 1.0e-14_dp .and. &
        maxval(abs(ledger - [0.64_dp, 1.47_dp, 1.65_dp, 0.675_dp])) &
        < 1.0e-14_dp .and. abs(total - 4.435_dp) < 1.0e-14_dp, &
        "FCI coupled split matches the independent additive work oracle")

    call evaluate_fci_coupled_field_split_ledger_jvp( &
        residual, parallel_action, plane_action, retained_actions, weights, &
        residual_dot, parallel_action_dot, plane_action_dot, &
        retained_actions_dot, weights_dot, correction_dot, ledger_dot, &
        total_dot, status)
    call evaluate_fci_coupled_field_split_ledger( &
        residual + step*residual_dot, &
        parallel_action + step*parallel_action_dot, &
        plane_action + step*plane_action_dot, &
        retained_actions + step*retained_actions_dot, weights + step*weights_dot, &
        correction_plus, ledger_plus, total_plus, status)
    call evaluate_fci_coupled_field_split_ledger( &
        residual - step*residual_dot, &
        parallel_action - step*parallel_action_dot, &
        plane_action - step*plane_action_dot, &
        retained_actions - step*retained_actions_dot, weights - step*weights_dot, &
        correction_minus, ledger_minus, total_minus, status)
    call check_condition( &
        maxval(abs(correction_dot - &
            (correction_plus - correction_minus)/(2.0_dp*step))) &
        < 2.0e-10_dp .and. &
        maxval(abs(ledger_dot - (ledger_plus - ledger_minus)/(2.0_dp*step))) &
        < 2.0e-10_dp .and. &
        abs(total_dot - (total_plus - total_minus)/(2.0_dp*step)) &
        < 2.0e-10_dp, &
        "FCI coupled split JVP matches a centered finite-difference oracle")

    call evaluate_fci_coupled_field_split_ledger_vjp( &
        residual, parallel_action, plane_action, retained_actions, weights, &
        correction_bar, ledger_bar, total_bar, residual_bar, &
        parallel_action_bar, plane_action_bar, retained_actions_bar, &
        weights_bar, status)
    lhs = dot_product(correction_bar, correction_dot) + &
        dot_product(ledger_bar, ledger_dot) + total_bar*total_dot
    rhs = dot_product(residual_bar, residual_dot) + &
        dot_product(parallel_action_bar, parallel_action_dot) + &
        dot_product(plane_action_bar, plane_action_dot) + &
        sum(retained_actions_bar*retained_actions_dot) + &
        dot_product(weights_bar, weights_dot)
    call check_condition(abs(lhs - rhs) < 2.0e-13_dp, &
        "FCI coupled split VJP satisfies the full adjoint identity")

    call evaluate_fci_coupled_field_split_ledger( &
        residual, -parallel_action, plane_action, retained_actions, weights, &
        correction, ledger, total, status)
    call check_condition(status%code /= 0 .and. &
        maxval(abs(correction)) == 0.0_dp .and. &
        maxval(abs(ledger)) == 0.0_dp .and. total == 0.0_dp, &
        "FCI coupled split rejects a negative-energy component")

    call evaluate_fci_coupled_field_split_ledger( &
        residual, parallel_action, plane_action, retained_actions, &
        [0.4_dp, 0.0_dp, 1.1_dp, 0.3_dp], correction, ledger, total, status)
    call check_condition(status%code /= 0, &
        "FCI coupled split rejects a non-positive partition weight")
    call check_summary("FCI coupled field-split ledger")
end program test_fci_coupled_field_split_ledger
