program test_surface_shape_objective
    !! Independent weighted fixed-topology surface-shape objective oracle.
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_surface_shape_objective, &
        evaluate_surface_shape_objective_jvp, &
        evaluate_surface_shape_objective_vjp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64
    real(dp), parameter :: finite_difference_step = 1.0e-7_dp
    real(dp) :: candidate(2, 3), target(2, 3), candidate_dot(2, 3)
    real(dp) :: target_dot(2, 3), weights(3), weights_dot(3)
    real(dp) :: candidate_plus(2, 3), candidate_minus(2, 3)
    real(dp) :: target_plus(2, 3), target_minus(2, 3)
    real(dp) :: weights_plus(3), weights_minus(3)
    real(dp) :: candidate_bar(2, 3), target_bar(2, 3), weights_bar(3)
    real(dp) :: objective, objective_dot, objective_plus, objective_minus
    real(dp) :: expected_objective, expected_dot, finite_difference_dot
    real(dp) :: adjoint_left
    real(dp) :: invalid_target(2, 2)
    type(fortsparse_status_t) :: status

    candidate = reshape([ &
        1.2_dp, -0.4_dp, 2.1_dp, 0.6_dp, 0.7_dp, -1.1_dp], [2, 3])
    target = reshape([ &
        1.0_dp, -0.5_dp, 2.0_dp, 0.5_dp, 0.9_dp, -1.0_dp], [2, 3])
    candidate_dot = reshape([ &
        0.3_dp, -0.2_dp, 0.1_dp, 0.4_dp, -0.1_dp, 0.2_dp], [2, 3])
    target_dot = reshape([ &
        -0.1_dp, 0.2_dp, 0.0_dp, -0.3_dp, 0.2_dp, -0.1_dp], [2, 3])
    weights = [1.0_dp, 2.0_dp, 0.5_dp]
    weights_dot = [0.2_dp, -0.1_dp, 0.3_dp]

    call evaluate_surface_shape_objective( &
        candidate, target, weights, objective, status)
    call check_condition(status%code == FORTSPARSE_OK, &
        "surface shape objective accepts finite fixed-topology samples")
    expected_objective = 0.5_dp*(weights(1)*((0.2_dp)**2 + 0.1_dp**2) + &
        weights(2)*((0.1_dp)**2 + 0.1_dp**2) + &
        weights(3)*((-0.2_dp)**2 + (-0.1_dp)**2))
    call check_condition(abs(objective - expected_objective) < 1.0e-14_dp, &
        "surface shape objective matches weighted Euclidean oracle")

    call evaluate_surface_shape_objective_jvp( &
        candidate, target, weights, candidate_dot, target_dot, weights_dot, &
        objective_dot, status)
    expected_dot = weights(1)*((0.2_dp)*(0.3_dp + 0.1_dp) + &
        0.1_dp*(-0.2_dp - 0.2_dp)) + &
        weights(2)*((0.1_dp)*(0.1_dp + 0.0_dp) + &
        0.1_dp*(0.4_dp + 0.3_dp)) + &
        weights(3)*((-0.2_dp)*(-0.1_dp - 0.2_dp) + &
        (-0.1_dp)*(0.2_dp + 0.1_dp)) + &
        0.5_dp*(weights_dot(1)*((0.2_dp)**2 + 0.1_dp**2) + &
        weights_dot(2)*((0.1_dp)**2 + 0.1_dp**2) + &
        weights_dot(3)*((-0.2_dp)**2 + (-0.1_dp)**2))
    call check_condition(status%code == FORTSPARSE_OK .and. &
        abs(objective_dot - expected_dot) < 1.0e-14_dp, &
        "surface shape JVP matches independent product rule")

    candidate_plus = candidate + finite_difference_step*candidate_dot
    candidate_minus = candidate - finite_difference_step*candidate_dot
    target_plus = target + finite_difference_step*target_dot
    target_minus = target - finite_difference_step*target_dot
    weights_plus = weights + finite_difference_step*weights_dot
    weights_minus = weights - finite_difference_step*weights_dot
    call evaluate_surface_shape_objective( &
        candidate_plus, target_plus, weights_plus, objective_plus, status)
    call evaluate_surface_shape_objective( &
        candidate_minus, target_minus, weights_minus, objective_minus, status)
    finite_difference_dot = (objective_plus - objective_minus)/ &
        (2.0_dp*finite_difference_step)
    call check_condition(abs(objective_dot - finite_difference_dot) < 2.0e-9_dp, &
        "surface shape JVP matches centered difference")

    call evaluate_surface_shape_objective_vjp( &
        candidate, target, weights, 1.0_dp, candidate_bar, target_bar, &
        weights_bar, status)
    adjoint_left = sum(candidate_bar*candidate_dot) + &
        sum(target_bar*target_dot) + dot_product(weights_bar, weights_dot)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        abs(adjoint_left - objective_dot) < 1.0e-14_dp, &
        "surface shape VJP satisfies real adjoint identity")

    weights(2) = -weights(2)
    call evaluate_surface_shape_objective( &
        candidate, target, weights, objective, status)
    call check_condition(status%code == FORTSPARSE_INVALID_MATRIX, &
        "surface shape objective rejects non-positive weights")
    weights(2) = 2.0_dp

    invalid_target = 0.0_dp
    call evaluate_surface_shape_objective( &
        candidate, invalid_target, weights, objective, status)
    call check_condition(status%code == FORTSPARSE_INVALID_MATRIX, &
        "surface shape objective rejects incompatible topology samples")

    call check_summary("surface shape objective")
end program test_surface_shape_objective
