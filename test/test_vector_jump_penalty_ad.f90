program test_vector_jump_penalty_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_vector_jump_penalty, &
        assemble_vector_jump_penalty_jvp, assemble_vector_jump_penalty_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: plus_trace(2, 2, 2) = reshape([ &
        1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 0.5_dp, -0.7_dp, 0.2_dp, 0.8_dp], [2, 2, 2])
    real(dp), parameter :: minus_trace(2, 1, 2) = reshape([ &
        5.0_dp, 6.0_dp, -0.4_dp, 0.9_dp], [2, 1, 2])
    real(dp), parameter :: surface_weights(2) = [2.0_dp, 1.0_dp]
    real(dp), parameter :: penalty = 2.5_dp
    real(dp), parameter :: component_metric(2, 2, 2) = reshape([ &
        2.0_dp, 0.3_dp, 0.3_dp, 1.5_dp, 1.2_dp, -0.2_dp, &
        0.4_dp, 0.8_dp], [2, 2, 2])
    real(dp), parameter :: plus_trace_dot(2, 2, 2) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp, -0.2_dp, 0.1_dp, 0.05_dp, 0.2_dp], [2, 2, 2])
    real(dp), parameter :: minus_trace_dot(2, 1, 2) = reshape([ &
        0.2_dp, -0.3_dp, 0.4_dp, -0.1_dp], [2, 1, 2])
    real(dp), parameter :: surface_weights_dot(2) = [0.05_dp, -0.1_dp]
    real(dp), parameter :: penalty_dot = 0.7_dp
    real(dp), parameter :: component_metric_dot(2, 2, 2) = reshape([ &
        -0.1_dp, 0.2_dp, 0.05_dp, -0.4_dp, 0.3_dp, 0.1_dp, &
        -0.2_dp, 0.15_dp], [2, 2, 2])
    real(dp), parameter :: matrix_bar(3, 3) = reshape([ &
        0.4_dp, -0.2_dp, 0.1_dp, 0.3_dp, 0.5_dp, -0.6_dp, &
        -0.7_dp, 0.8_dp, 0.9_dp], [3, 3])
    real(dp) :: matrix(3, 3), matrix_dot(3, 3), matrix_oracle(3, 3)
    real(dp) :: plus_trace_bar(2, 2, 2), minus_trace_bar(2, 1, 2)
    real(dp) :: surface_weights_bar(2), penalty_bar
    real(dp) :: component_metric_bar(2, 2, 2), lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_vector_jump_penalty( &
        plus_trace, minus_trace, surface_weights, penalty, component_metric, &
        matrix, status)
    call oracle_vector_penalty(matrix_oracle)
    call check_condition(status%code == 0 .and. &
        maxval(abs(matrix - matrix_oracle)) < 1.0e-14_dp, &
        "vector jump penalty matches the independent tensor oracle")

    call assemble_vector_jump_penalty_jvp( &
        plus_trace, minus_trace, surface_weights, penalty, component_metric, &
        plus_trace_dot, minus_trace_dot, surface_weights_dot, penalty_dot, &
        component_metric_dot, matrix_dot, status)
    call check_condition(status%code == 0, &
        "vector jump penalty JVP accepts trace, metric, and weight directions")

    call assemble_vector_jump_penalty_vjp( &
        plus_trace, minus_trace, surface_weights, penalty, component_metric, &
        matrix_bar, plus_trace_bar, minus_trace_bar, surface_weights_bar, &
        penalty_bar, component_metric_bar, status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(plus_trace_bar*plus_trace_dot) + sum(minus_trace_bar*minus_trace_dot) + &
        sum(surface_weights_bar*surface_weights_dot) + penalty_bar*penalty_dot + &
        sum(component_metric_bar*component_metric_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "vector jump penalty VJP satisfies the real dot-product identity")

    call assemble_vector_jump_penalty_jvp( &
        plus_trace, minus_trace, surface_weights, penalty, component_metric, &
        plus_trace(:, 1:1, :), minus_trace_dot, surface_weights_dot, penalty_dot, &
        component_metric_dot, matrix_dot, status)
    call check_condition(status%code /= 0, &
        "vector jump penalty JVP rejects incompatible trace dimensions")
    call check_summary("vector jump penalty AD")

contains

    subroutine oracle_vector_penalty(result)
        real(dp), intent(out) :: result(:, :)
        integer :: q, row, column, a, b
        real(dp) :: jump(3, 2)

        result = 0.0_dp
        do q = 1, 2
            jump(1, :) = plus_trace(q, 1, :)
            jump(2, :) = plus_trace(q, 2, :)
            jump(3, :) = -minus_trace(q, 1, :)
            do row = 1, 3
                do column = 1, 3
                    do a = 1, 2
                        do b = 1, 2
                            result(row, column) = result(row, column) + &
                                penalty*surface_weights(q)*jump(row, a)* &
                                component_metric(q, a, b)*jump(column, b)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine oracle_vector_penalty

end program test_vector_jump_penalty_ad
