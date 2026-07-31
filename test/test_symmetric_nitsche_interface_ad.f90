program test_symmetric_nitsche_interface_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_symmetric_nitsche_interface, &
        assemble_symmetric_nitsche_interface_jvp, &
        assemble_symmetric_nitsche_interface_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: plus_trace(2, 2) = reshape([ &
        1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2, 2])
    real(dp), parameter :: minus_trace(2, 1) = reshape([5.0_dp, 6.0_dp], [2, 1])
    real(dp), parameter :: plus_flux(2, 2) = reshape([ &
        0.5_dp, -0.2_dp, 0.7_dp, 0.4_dp], [2, 2])
    real(dp), parameter :: minus_flux(2, 1) = reshape([0.3_dp, -0.6_dp], [2, 1])
    real(dp), parameter :: surface_weights(2) = [2.0_dp, 1.0_dp]
    real(dp), parameter :: penalty = 4.0_dp
    real(dp), parameter :: plus_trace_dot(2, 2) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp], [2, 2])
    real(dp), parameter :: minus_trace_dot(2, 1) = reshape([0.2_dp, -0.3_dp], [2, 1])
    real(dp), parameter :: plus_flux_dot(2, 2) = reshape([ &
        -0.1_dp, 0.2_dp, 0.05_dp, -0.4_dp], [2, 2])
    real(dp), parameter :: minus_flux_dot(2, 1) = reshape([0.3_dp, 0.1_dp], [2, 1])
    real(dp), parameter :: surface_weights_dot(2) = [0.05_dp, -0.1_dp]
    real(dp), parameter :: penalty_dot = 0.7_dp
    real(dp), parameter :: matrix_bar(3, 3) = reshape([ &
        0.4_dp, -0.2_dp, 0.1_dp, 0.3_dp, 0.5_dp, -0.6_dp, &
        -0.7_dp, 0.8_dp, 0.9_dp], [3, 3])
    real(dp) :: matrix(3, 3), matrix_dot(3, 3)
    real(dp) :: plus_trace_bar(2, 2), minus_trace_bar(2, 1)
    real(dp) :: plus_flux_bar(2, 2), minus_flux_bar(2, 1)
    real(dp) :: surface_weights_bar(2), penalty_bar, lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_symmetric_nitsche_interface( &
        plus_trace, minus_trace, plus_flux, minus_flux, surface_weights, penalty, &
        matrix, status)
    call check_condition(status%code == 0, &
        "Nitsche AD test accepts the primal block")
    call assemble_symmetric_nitsche_interface_jvp( &
        plus_trace, minus_trace, plus_flux, minus_flux, surface_weights, penalty, &
        plus_trace_dot, minus_trace_dot, plus_flux_dot, minus_flux_dot, &
        surface_weights_dot, penalty_dot, matrix_dot, status)
    call check_condition(status%code == 0, &
        "Nitsche JVP accepts trace, flux, weight, and penalty directions")

    call assemble_symmetric_nitsche_interface_vjp( &
        plus_trace, minus_trace, plus_flux, minus_flux, surface_weights, penalty, &
        matrix_bar, plus_trace_bar, minus_trace_bar, plus_flux_bar, minus_flux_bar, &
        surface_weights_bar, penalty_bar, status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(plus_trace_bar*plus_trace_dot) + sum(minus_trace_bar*minus_trace_dot) + &
        sum(plus_flux_bar*plus_flux_dot) + sum(minus_flux_bar*minus_flux_dot) + &
        sum(surface_weights_bar*surface_weights_dot) + penalty_bar*penalty_dot
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "Nitsche VJP satisfies the real dot-product identity")

    call assemble_symmetric_nitsche_interface_jvp( &
        plus_trace, minus_trace, plus_flux, minus_flux, surface_weights, penalty, &
        plus_trace_dot(:, 1:1), minus_trace_dot, plus_flux_dot, minus_flux_dot, &
        surface_weights_dot, penalty_dot, matrix_dot, status)
    call check_condition(status%code /= 0, &
        "Nitsche JVP rejects incompatible trace dimensions")
    call check_summary("symmetric Nitsche interface AD")
end program test_symmetric_nitsche_interface_ad
