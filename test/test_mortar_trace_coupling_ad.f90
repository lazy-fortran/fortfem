program test_mortar_trace_coupling_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_mortar_trace_coupling, &
        assemble_mortar_trace_coupling_jvp, assemble_mortar_trace_coupling_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: test_trace(2, 2) = reshape([ &
        1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2, 2])
    real(dp), parameter :: trial_trace(2, 1) = reshape([5.0_dp, 6.0_dp], [2, 1])
    real(dp), parameter :: surface_weights(2) = [2.0_dp, 1.0_dp]
    real(dp), parameter :: test_trace_dot(2, 2) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp], [2, 2])
    real(dp), parameter :: trial_trace_dot(2, 1) = reshape([0.2_dp, -0.3_dp], [2, 1])
    real(dp), parameter :: surface_weights_dot(2) = [0.05_dp, -0.1_dp]
    real(dp), parameter :: matrix_bar(2, 1) = reshape([0.7_dp, -0.4_dp], [2, 1])
    real(dp) :: matrix(2, 1), matrix_dot(2, 1)
    real(dp) :: test_trace_bar(2, 2), trial_trace_bar(2, 1)
    real(dp) :: surface_weights_bar(2), lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_mortar_trace_coupling( &
        test_trace, trial_trace, surface_weights, matrix, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(matrix - reshape([22.0_dp, 54.0_dp], [2, 1]))) < 1.0e-14_dp, &
        "mortar AD test starts from the independent cross-mass oracle")

    call assemble_mortar_trace_coupling_jvp( &
        test_trace, trial_trace, surface_weights, test_trace_dot, trial_trace_dot, &
        surface_weights_dot, matrix_dot, status)
    call check_condition(status%code == 0, &
        "mortar cross-mass JVP accepts compatible directions")

    call assemble_mortar_trace_coupling_vjp( &
        test_trace, trial_trace, surface_weights, matrix_bar, test_trace_bar, &
        trial_trace_bar, surface_weights_bar, status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(test_trace_bar*test_trace_dot) + sum(trial_trace_bar*trial_trace_dot) + &
        sum(surface_weights_bar*surface_weights_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-14_dp, &
        "mortar cross-mass VJP satisfies the real dot-product identity")

    call assemble_mortar_trace_coupling_jvp( &
        test_trace, trial_trace, surface_weights, test_trace_dot(:, 1:1), &
        trial_trace_dot, surface_weights_dot, matrix_dot, status)
    call check_condition(status%code /= 0, &
        "mortar cross-mass JVP rejects incompatible trace dimensions")
    call check_summary("mortar trace coupling AD")
end program test_mortar_trace_coupling_ad
