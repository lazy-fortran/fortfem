program test_interface_jump_penalty_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_interface_jump_penalty, &
        assemble_interface_jump_penalty_jvp, assemble_interface_jump_penalty_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: plus_trace(2, 2) = reshape([ &
        1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2, 2])
    real(dp), parameter :: minus_trace(2, 1) = reshape([5.0_dp, 6.0_dp], [2, 1])
    real(dp), parameter :: surface_weights(2) = [2.0_dp, 1.0_dp]
    real(dp), parameter :: penalty = 4.0_dp
    real(dp), parameter :: plus_trace_dot(2, 2) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp], [2, 2])
    real(dp), parameter :: minus_trace_dot(2, 1) = reshape([0.2_dp, -0.3_dp], [2, 1])
    real(dp), parameter :: surface_weights_dot(2) = [0.05_dp, -0.1_dp]
    real(dp), parameter :: penalty_dot = 0.7_dp
    real(dp), parameter :: matrix_bar(3, 3) = reshape([ &
        0.4_dp, -0.2_dp, 0.1_dp, 0.3_dp, 0.5_dp, -0.6_dp, &
        -0.7_dp, 0.8_dp, 0.9_dp], [3, 3])
    real(dp) :: matrix(3, 3), matrix_dot(3, 3)
    real(dp) :: plus_trace_bar(2, 2), minus_trace_bar(2, 1)
    real(dp) :: surface_weights_bar(2), penalty_bar, lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_interface_jump_penalty( &
        plus_trace, minus_trace, surface_weights, penalty, matrix, status)
    call check_condition(status%code == 0, &
        "jump penalty AD test accepts the primal block")
    call assemble_interface_jump_penalty_jvp( &
        plus_trace, minus_trace, surface_weights, penalty, plus_trace_dot, &
        minus_trace_dot, surface_weights_dot, penalty_dot, matrix_dot, status)
    call check_condition(status%code == 0, &
        "jump penalty JVP accepts trace, weight, and penalty directions")

    call assemble_interface_jump_penalty_vjp( &
        plus_trace, minus_trace, surface_weights, penalty, matrix_bar, &
        plus_trace_bar, minus_trace_bar, surface_weights_bar, penalty_bar, status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(plus_trace_bar*plus_trace_dot) + sum(minus_trace_bar*minus_trace_dot) + &
        sum(surface_weights_bar*surface_weights_dot) + penalty_bar*penalty_dot
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "jump penalty VJP satisfies the real dot-product identity")

    call assemble_interface_jump_penalty_jvp( &
        plus_trace, minus_trace, surface_weights, penalty, plus_trace_dot(:, 1:1), &
        minus_trace_dot, surface_weights_dot, penalty_dot, matrix_dot, status)
    call check_condition(status%code /= 0, &
        "jump penalty JVP rejects incompatible trace dimensions")
    call check_summary("interface jump penalty AD")
end program test_interface_jump_penalty_ad
