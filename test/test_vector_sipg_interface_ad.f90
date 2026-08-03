program test_vector_sipg_interface_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_vector_sipg_interface, &
        assemble_vector_sipg_interface_jvp, assemble_vector_sipg_interface_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: test_plus_trace(2, 2, 2) = reshape([ &
        1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 0.5_dp, -0.7_dp, 0.2_dp, 0.8_dp], [2, 2, 2])
    real(dp), parameter :: test_minus_trace(2, 1, 2) = reshape([ &
        5.0_dp, 6.0_dp, -0.4_dp, 0.9_dp], [2, 1, 2])
    real(dp), parameter :: trial_plus_trace(2, 1, 2) = reshape([ &
        0.3_dp, -0.1_dp, 0.4_dp, 0.8_dp], [2, 1, 2])
    real(dp), parameter :: trial_minus_trace(2, 2, 2) = reshape([ &
        0.2_dp, 0.7_dp, -0.4_dp, 0.5_dp, 0.6_dp, -0.2_dp, &
        0.1_dp, 0.9_dp], [2, 2, 2])
    real(dp), parameter :: test_plus_flux(2, 2, 2) = reshape([ &
        0.5_dp, -0.2_dp, 0.7_dp, 0.4_dp, 0.2_dp, 0.9_dp, &
        -0.1_dp, 0.6_dp], [2, 2, 2])
    real(dp), parameter :: test_minus_flux(2, 1, 2) = reshape([ &
        0.3_dp, -0.6_dp, -0.2_dp, 0.4_dp], [2, 1, 2])
    real(dp), parameter :: trial_plus_flux(2, 1, 2) = reshape([ &
        0.2_dp, 0.9_dp, 0.1_dp, -0.3_dp], [2, 1, 2])
    real(dp), parameter :: trial_minus_flux(2, 2, 2) = reshape([ &
        -0.1_dp, 0.6_dp, 0.4_dp, -0.3_dp, 0.7_dp, 0.2_dp, &
        -0.5_dp, 0.8_dp], [2, 2, 2])
    real(dp), parameter :: surface_weights(2) = [2.0_dp, 1.0_dp]
    real(dp), parameter :: penalty = 2.5_dp
    integer, parameter :: consistency_sign = 1
    real(dp), parameter :: component_metric(2, 2, 2) = reshape([ &
        2.0_dp, 0.3_dp, 0.3_dp, 1.5_dp, 1.2_dp, -0.2_dp, &
        0.4_dp, 0.8_dp], [2, 2, 2])
    real(dp), parameter :: test_plus_trace_dot(2, 2, 2) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp, -0.2_dp, 0.1_dp, &
        0.05_dp, 0.2_dp], [2, 2, 2])
    real(dp), parameter :: test_minus_trace_dot(2, 1, 2) = reshape([ &
        0.2_dp, -0.3_dp, 0.4_dp, -0.1_dp], [2, 1, 2])
    real(dp), parameter :: trial_plus_trace_dot(2, 1, 2) = reshape([ &
        0.05_dp, 0.2_dp, -0.1_dp, 0.3_dp], [2, 1, 2])
    real(dp), parameter :: trial_minus_trace_dot(2, 2, 2) = reshape([ &
        -0.1_dp, 0.3_dp, 0.4_dp, 0.2_dp, 0.1_dp, -0.2_dp, &
        0.3_dp, 0.4_dp], [2, 2, 2])
    real(dp), parameter :: test_plus_flux_dot(2, 2, 2) = reshape([ &
        -0.1_dp, 0.2_dp, 0.05_dp, -0.4_dp, 0.3_dp, 0.1_dp, &
        -0.2_dp, 0.15_dp], [2, 2, 2])
    real(dp), parameter :: test_minus_flux_dot(2, 1, 2) = reshape([ &
        0.3_dp, 0.1_dp, -0.2_dp, 0.25_dp], [2, 1, 2])
    real(dp), parameter :: trial_plus_flux_dot(2, 1, 2) = reshape([ &
        0.2_dp, -0.15_dp, 0.1_dp, 0.4_dp], [2, 1, 2])
    real(dp), parameter :: trial_minus_flux_dot(2, 2, 2) = reshape([ &
        0.1_dp, -0.3_dp, 0.2_dp, 0.4_dp, -0.2_dp, 0.1_dp, &
        0.3_dp, -0.1_dp], [2, 2, 2])
    real(dp), parameter :: surface_weights_dot(2) = [0.05_dp, -0.1_dp]
    real(dp), parameter :: penalty_dot = 0.7_dp
    real(dp), parameter :: component_metric_dot(2, 2, 2) = reshape([ &
        -0.1_dp, 0.2_dp, 0.05_dp, -0.4_dp, 0.3_dp, 0.1_dp, &
        -0.2_dp, 0.15_dp], [2, 2, 2])
    real(dp), parameter :: matrix_bar(3, 3) = reshape([ &
        0.4_dp, -0.2_dp, 0.1_dp, 0.3_dp, 0.5_dp, -0.6_dp, &
        -0.7_dp, 0.8_dp, 0.9_dp], [3, 3])
    real(dp) :: matrix(3, 3), matrix_dot(3, 3), matrix_oracle(3, 3)
    real(dp) :: test_plus_trace_bar(2, 2, 2), test_minus_trace_bar(2, 1, 2)
    real(dp) :: trial_plus_trace_bar(2, 1, 2), trial_minus_trace_bar(2, 2, 2)
    real(dp) :: test_plus_flux_bar(2, 2, 2), test_minus_flux_bar(2, 1, 2)
    real(dp) :: trial_plus_flux_bar(2, 1, 2), trial_minus_flux_bar(2, 2, 2)
    real(dp) :: surface_weights_bar(2), penalty_bar, component_metric_bar(2, 2, 2)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_vector_sipg_interface( &
        test_plus_trace, test_minus_trace, trial_plus_trace, trial_minus_trace, &
        test_plus_flux, test_minus_flux, trial_plus_flux, trial_minus_flux, &
        surface_weights, penalty, consistency_sign, component_metric, matrix, status)
    call oracle_vector_sipg(matrix_oracle)
    call check_condition(status%code == 0 .and. &
        maxval(abs(matrix - matrix_oracle)) < 1.0e-14_dp, &
        "vector SIPG block matches the independent metric oracle")

    call assemble_vector_sipg_interface_jvp( &
        test_plus_trace, test_minus_trace, trial_plus_trace, trial_minus_trace, &
        test_plus_flux, test_minus_flux, trial_plus_flux, trial_minus_flux, &
        surface_weights, penalty, consistency_sign, component_metric, &
        test_plus_trace_dot, test_minus_trace_dot, trial_plus_trace_dot, &
        trial_minus_trace_dot, test_plus_flux_dot, test_minus_flux_dot, &
        trial_plus_flux_dot, trial_minus_flux_dot, surface_weights_dot, penalty_dot, &
        component_metric_dot, matrix_dot, status)
    call check_condition(status%code == 0, &
        "vector SIPG JVP accepts all trace, flux, metric, and weight directions")

    call assemble_vector_sipg_interface_vjp( &
        test_plus_trace, test_minus_trace, trial_plus_trace, trial_minus_trace, &
        test_plus_flux, test_minus_flux, trial_plus_flux, trial_minus_flux, &
        surface_weights, penalty, consistency_sign, component_metric, matrix_bar, &
        test_plus_trace_bar, test_minus_trace_bar, trial_plus_trace_bar, &
        trial_minus_trace_bar, test_plus_flux_bar, test_minus_flux_bar, &
        trial_plus_flux_bar, trial_minus_flux_bar, surface_weights_bar, penalty_bar, &
        component_metric_bar, status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(test_plus_trace_bar*test_plus_trace_dot) + &
        sum(test_minus_trace_bar*test_minus_trace_dot) + &
        sum(trial_plus_trace_bar*trial_plus_trace_dot) + &
        sum(trial_minus_trace_bar*trial_minus_trace_dot) + &
        sum(test_plus_flux_bar*test_plus_flux_dot) + &
        sum(test_minus_flux_bar*test_minus_flux_dot) + &
        sum(trial_plus_flux_bar*trial_plus_flux_dot) + &
        sum(trial_minus_flux_bar*trial_minus_flux_dot) + &
        sum(surface_weights_bar*surface_weights_dot) + penalty_bar*penalty_dot + &
        sum(component_metric_bar*component_metric_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "vector SIPG VJP satisfies the real dot-product identity")

    call assemble_vector_sipg_interface( &
        test_plus_trace, test_minus_trace, trial_plus_trace, trial_minus_trace, &
        test_plus_flux, test_minus_flux, trial_plus_flux, trial_minus_flux, &
        surface_weights, penalty, 2, component_metric, matrix, status)
    call check_condition(status%code /= 0, &
        "vector SIPG rejects a consistency selector outside the standard variants")
    call check_summary("vector SIPG interface AD")

contains

    subroutine oracle_vector_sipg(result)
        real(dp), intent(out) :: result(:, :)
        integer :: q, row, column, a, b
        real(dp) :: test_jump(3, 2), trial_jump(3, 2)
        real(dp) :: test_average(3, 2), trial_average(3, 2)

        result = 0.0_dp
        do q = 1, 2
            test_jump(1, :) = test_plus_trace(q, 1, :)
            test_jump(2, :) = test_plus_trace(q, 2, :)
            test_jump(3, :) = -test_minus_trace(q, 1, :)
            trial_jump(1, :) = trial_plus_trace(q, 1, :)
            trial_jump(2, :) = -trial_minus_trace(q, 1, :)
            trial_jump(3, :) = -trial_minus_trace(q, 2, :)
            test_average(1, :) = 0.5_dp*test_plus_flux(q, 1, :)
            test_average(2, :) = 0.5_dp*test_plus_flux(q, 2, :)
            test_average(3, :) = 0.5_dp*test_minus_flux(q, 1, :)
            trial_average(1, :) = 0.5_dp*trial_plus_flux(q, 1, :)
            trial_average(2, :) = 0.5_dp*trial_minus_flux(q, 1, :)
            trial_average(3, :) = 0.5_dp*trial_minus_flux(q, 2, :)
            do row = 1, 3
                do column = 1, 3
                    do a = 1, 2
                        do b = 1, 2
                            result(row, column) = result(row, column) + &
                                surface_weights(q)*( &
                                -test_average(row, a)*component_metric(q, a, b)* &
                                trial_jump(column, b) - consistency_sign* &
                                test_jump(row, a)*component_metric(q, a, b)* &
                                trial_average(column, b) + penalty*test_jump(row, a)* &
                                component_metric(q, a, b)*trial_jump(column, b))
                        end do
                    end do
                end do
            end do
        end do
    end subroutine oracle_vector_sipg

end program test_vector_sipg_interface_ad
