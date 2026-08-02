program test_surface_integral_constraint
    !! Independent fixed-topology weighted surface-integral constraint oracle.
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_surface_integral_constraint, &
        evaluate_surface_integral_constraint_jvp, &
        evaluate_surface_integral_constraint_vjp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64
    real(dp), parameter :: finite_difference_step = 1.0e-7_dp
    real(dp) :: samples(3), weights(3), samples_dot(3), weights_dot(3)
    real(dp) :: samples_plus(3), samples_minus(3), weights_plus(3)
    real(dp) :: weights_minus(3), samples_bar(3), weights_bar(3)
    real(dp) :: target, target_dot, target_plus, target_minus, target_bar
    real(dp) :: constraint, constraint_dot, constraint_plus, constraint_minus
    real(dp) :: expected_constraint, expected_dot, finite_difference_dot
    real(dp) :: constraint_bar, adjoint_left
    type(fortsparse_status_t) :: status

    samples = [2.0_dp, -1.0_dp, 0.5_dp]
    weights = [1.0_dp, 2.0_dp, 0.5_dp]
    target = 0.75_dp
    samples_dot = [0.3_dp, -0.2_dp, 0.4_dp]
    weights_dot = [0.1_dp, -0.1_dp, 0.2_dp]
    target_dot = 0.05_dp

    call evaluate_surface_integral_constraint( &
        samples, weights, target, constraint, status)
    call check_condition(status%code == FORTSPARSE_OK, &
        "surface-integral constraint accepts finite samples and weights")
    expected_constraint = (1.0_dp*2.0_dp + 2.0_dp*(-1.0_dp) + &
        0.5_dp*0.5_dp) - 0.75_dp
    call check_condition(abs(constraint - expected_constraint) < 1.0e-14_dp, &
        "surface-integral constraint matches independent weighted ledger")

    call evaluate_surface_integral_constraint_jvp( &
        samples, weights, target, samples_dot, weights_dot, target_dot, &
        constraint_dot, status)
    expected_dot = 0.1_dp*2.0_dp + 1.0_dp*0.3_dp + &
        (-0.1_dp)*(-1.0_dp) + 2.0_dp*(-0.2_dp) + &
        0.2_dp*0.5_dp + 0.5_dp*0.4_dp - 0.05_dp
    call check_condition(status%code == FORTSPARSE_OK .and. &
        abs(constraint_dot - expected_dot) < 1.0e-14_dp, &
        "surface-integral JVP matches independent product rule")

    samples_plus = samples + finite_difference_step*samples_dot
    samples_minus = samples - finite_difference_step*samples_dot
    weights_plus = weights + finite_difference_step*weights_dot
    weights_minus = weights - finite_difference_step*weights_dot
    target_plus = target + finite_difference_step*target_dot
    target_minus = target - finite_difference_step*target_dot
    call evaluate_surface_integral_constraint( &
        samples_plus, weights_plus, target_plus, constraint_plus, status)
    call evaluate_surface_integral_constraint( &
        samples_minus, weights_minus, target_minus, constraint_minus, status)
    finite_difference_dot = (constraint_plus - constraint_minus)/ &
        (2.0_dp*finite_difference_step)
    call check_condition(abs(constraint_dot - finite_difference_dot) < 1.0e-8_dp, &
        "surface-integral JVP matches centered difference")

    constraint_bar = 1.3_dp
    call evaluate_surface_integral_constraint_vjp( &
        samples, weights, target, constraint_bar, samples_bar, weights_bar, &
        target_bar, status)
    adjoint_left = dot_product(samples_bar, samples_dot) + &
        dot_product(weights_bar, weights_dot) + target_bar*target_dot
    call check_condition(status%code == FORTSPARSE_OK .and. &
        abs(adjoint_left - constraint_bar*constraint_dot) < 1.0e-14_dp, &
        "surface-integral VJP satisfies real adjoint identity")

    weights(2) = 0.0_dp
    call evaluate_surface_integral_constraint( &
        samples, weights, target, constraint, status)
    call check_condition(status%code == FORTSPARSE_INVALID_MATRIX, &
        "surface-integral constraint rejects non-positive weights")
    weights(2) = 2.0_dp

    call evaluate_surface_integral_constraint( &
        samples(1:2), weights, target, constraint, status)
    call check_condition(status%code == FORTSPARSE_INVALID_MATRIX, &
        "surface-integral constraint rejects incompatible topology samples")

    call check_summary("surface integral constraint")
end program test_surface_integral_constraint
