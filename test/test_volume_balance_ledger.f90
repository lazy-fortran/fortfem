program test_volume_balance_ledger
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_volume_balance_ledger, &
        assemble_volume_balance_ledger_jvp, assemble_volume_balance_ledger_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: cell_count = 3, component_count = 2
    real(dp), parameter :: cell_weights(cell_count) = [2.0_dp, 1.5_dp, 0.75_dp]
    real(dp), parameter :: source_rate(cell_count, component_count) = reshape([ &
        1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp, -1.0_dp, 4.0_dp], &
        [cell_count, component_count])
    real(dp), parameter :: cell_weights_dot(cell_count) = [0.1_dp, -0.2_dp, 0.05_dp]
    real(dp), parameter :: source_rate_dot(cell_count, component_count) = reshape([ &
        -0.2_dp, 0.4_dp, 0.3_dp, -0.1_dp, 0.5_dp, 0.2_dp], &
        [cell_count, component_count])
    real(dp), parameter :: cell_ledger_bar(cell_count, component_count) = reshape([ &
        0.2_dp, -0.4_dp, 0.6_dp, 0.8_dp, -0.3_dp, 0.5_dp], &
        [cell_count, component_count])
    real(dp), parameter :: global_ledger_bar(component_count) = [0.7_dp, -0.5_dp]
    real(dp), parameter :: step = 1.0e-7_dp

    real(dp) :: cell_ledger(cell_count, component_count)
    real(dp) :: global_ledger(component_count)
    real(dp) :: cell_ledger_expected(cell_count, component_count)
    real(dp) :: global_expected(component_count)
    real(dp) :: cell_ledger_dot(cell_count, component_count)
    real(dp) :: global_dot(component_count)
    real(dp) :: cell_ledger_plus(cell_count, component_count)
    real(dp) :: cell_ledger_minus(cell_count, component_count)
    real(dp) :: global_plus(component_count), global_minus(component_count)
    real(dp) :: cell_weights_bar(cell_count)
    real(dp) :: source_rate_bar(cell_count, component_count)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_volume_balance_ledger( &
        cell_weights, source_rate, cell_ledger, global_ledger, status)
    call independent_value_oracle( &
        cell_weights, source_rate, cell_ledger_expected, global_expected)
    call check_condition(status%code == 0, &
        "volume balance ledger accepts positive cell measures")
    call check_condition(maxval(abs(cell_ledger - cell_ledger_expected)) < 1.0e-14_dp .and. &
        maxval(abs(global_ledger - global_expected)) < 1.0e-14_dp, &
        "volume balance ledger matches the independent integration oracle")

    call assemble_volume_balance_ledger_jvp( &
        cell_weights, source_rate, cell_weights_dot, source_rate_dot, &
        cell_ledger_dot, global_dot, status)
    call assemble_volume_balance_ledger( &
        cell_weights + step*cell_weights_dot, &
        source_rate + step*source_rate_dot, cell_ledger_plus, global_plus, status)
    call assemble_volume_balance_ledger( &
        cell_weights - step*cell_weights_dot, &
        source_rate - step*source_rate_dot, cell_ledger_minus, global_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(cell_ledger_dot - (cell_ledger_plus - cell_ledger_minus) / &
        (2.0_dp*step))) < 5.0e-9_dp .and. &
        maxval(abs(global_dot - (global_plus - global_minus) / &
        (2.0_dp*step))) < 5.0e-9_dp, &
        "volume balance ledger JVP matches central differences")

    call assemble_volume_balance_ledger_vjp( &
        cell_weights, source_rate, cell_ledger_bar, global_ledger_bar, &
        cell_weights_bar, source_rate_bar, status)
    lhs = sum(cell_ledger_bar*cell_ledger_dot) + sum(global_ledger_bar*global_dot)
    rhs = sum(cell_weights_bar*cell_weights_dot) + &
        sum(source_rate_bar*source_rate_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "volume balance ledger VJP satisfies the dot-product identity")

    call assemble_volume_balance_ledger( &
        [2.0_dp, 0.0_dp, 0.75_dp], source_rate, cell_ledger, global_ledger, status)
    call check_condition(status%code /= 0, &
        "volume balance ledger rejects non-positive cell measures")
    call check_summary("volume balance ledger")

contains

    subroutine independent_value_oracle(cell_weights, source_rate, cell_ledger, &
            global_ledger)
        real(dp), intent(in) :: cell_weights(:), source_rate(:, :)
        real(dp), intent(out) :: cell_ledger(:, :), global_ledger(:)
        integer :: cell, component

        cell_ledger = 0.0_dp
        global_ledger = 0.0_dp
        do cell = 1, size(cell_weights)
            do component = 1, size(source_rate, 2)
                cell_ledger(cell, component) = &
                    cell_weights(cell)*source_rate(cell, component)
                global_ledger(component) = global_ledger(component) + &
                    cell_ledger(cell, component)
            end do
        end do
    end subroutine independent_value_oracle

end program test_volume_balance_ledger
