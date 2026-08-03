program test_fci_terminal_boundary_ledger
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_fci_terminal_boundary_ledger, &
        assemble_fci_terminal_boundary_ledger_jvp, &
        assemble_fci_terminal_boundary_ledger_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: event_count = 3, cell_count = 2, component_count = 2
    integer, parameter :: owners(event_count) = [1, 2, 1]
    real(dp), parameter :: weights(event_count) = [2.0_dp, 1.5_dp, 0.5_dp]
    real(dp), parameter :: flux(event_count, component_count) = reshape([ &
        1.0_dp, -2.0_dp, 3.0_dp, 0.5_dp, -1.0_dp, 4.0_dp], &
        [event_count, component_count])
    real(dp), parameter :: volumes(cell_count) = [2.0_dp, 4.0_dp]
    real(dp), parameter :: weights_dot(event_count) = [0.1_dp, -0.2_dp, 0.3_dp]
    real(dp), parameter :: flux_dot(event_count, component_count) = reshape([ &
        -0.2_dp, 0.4_dp, 0.3_dp, -0.1_dp, 0.5_dp, 0.2_dp], &
        [event_count, component_count])
    real(dp), parameter :: volumes_dot(cell_count) = [0.2_dp, -0.3_dp]
    real(dp), parameter :: contribution_bar(cell_count, component_count) = reshape([ &
        0.2_dp, -0.4_dp, 0.6_dp, 0.8_dp], [cell_count, component_count])
    real(dp), parameter :: global_bar(component_count) = [0.7_dp, -0.5_dp]
    real(dp), parameter :: step = 1.0e-7_dp

    real(dp) :: contribution(cell_count, component_count)
    real(dp) :: global_ledger(component_count)
    real(dp) :: contribution_expected(cell_count, component_count)
    real(dp) :: global_expected(component_count)
    real(dp) :: contribution_dot(cell_count, component_count)
    real(dp) :: global_dot(component_count)
    real(dp) :: contribution_plus(cell_count, component_count)
    real(dp) :: contribution_minus(cell_count, component_count)
    real(dp) :: global_plus(component_count), global_minus(component_count)
    real(dp) :: weights_bar(event_count)
    real(dp) :: flux_bar(event_count, component_count)
    real(dp) :: volumes_bar(cell_count)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_fci_terminal_boundary_ledger( &
        owners, weights, flux, volumes, contribution, global_ledger, status)
    call independent_value_oracle(owners, weights, flux, volumes, &
        contribution_expected, global_expected)
    call check_condition(status%code == 0, &
        "vector FCI terminal ledger accepts positive event data")
    call check_condition(maxval(abs(contribution - contribution_expected)) < 1.0e-14_dp .and. &
        maxval(abs(global_ledger - global_expected)) < 1.0e-14_dp, &
        "vector FCI terminal ledger matches the independent tally oracle")

    call assemble_fci_terminal_boundary_ledger_jvp( &
        owners, weights, flux, volumes, weights_dot, flux_dot, volumes_dot, &
        contribution_dot, global_dot, status)
    call assemble_fci_terminal_boundary_ledger( &
        owners, weights + step*weights_dot, flux + step*flux_dot, &
        volumes + step*volumes_dot, contribution_plus, global_plus, status)
    call assemble_fci_terminal_boundary_ledger( &
        owners, weights - step*weights_dot, flux - step*flux_dot, &
        volumes - step*volumes_dot, contribution_minus, global_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(contribution_dot - (contribution_plus - contribution_minus) / &
        (2.0_dp*step))) < 5.0e-9_dp .and. &
        maxval(abs(global_dot - (global_plus - global_minus) / &
        (2.0_dp*step))) < 5.0e-9_dp, &
        "vector FCI terminal ledger JVP matches central differences")

    call assemble_fci_terminal_boundary_ledger_vjp( &
        owners, weights, flux, volumes, contribution_bar, global_bar, &
        weights_bar, flux_bar, volumes_bar, status)
    lhs = sum(contribution_bar*contribution_dot) + sum(global_bar*global_dot)
    rhs = sum(weights_bar*weights_dot) + sum(flux_bar*flux_dot) + &
        sum(volumes_bar*volumes_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "vector FCI terminal ledger VJP satisfies the dot-product identity")

    call assemble_fci_terminal_boundary_ledger( &
        owners, weights, flux, [2.0_dp, 0.0_dp], contribution, global_ledger, status)
    call check_condition(status%code /= 0, &
        "vector FCI terminal ledger rejects non-positive canonical volumes")
    call check_summary("vector FCI terminal ledger")

contains

    subroutine independent_value_oracle(owners, weights, flux, volumes, &
            contribution, global_ledger)
        integer, intent(in) :: owners(:)
        real(dp), intent(in) :: weights(:), flux(:, :), volumes(:)
        real(dp), intent(out) :: contribution(:, :), global_ledger(:)
        integer :: event, component, owner

        contribution = 0.0_dp
        global_ledger = 0.0_dp
        do event = 1, size(owners)
            owner = owners(event)
            do component = 1, size(flux, 2)
                contribution(owner, component) = contribution(owner, component) + &
                    weights(event)*flux(event, component)/volumes(owner)
                global_ledger(component) = global_ledger(component) + &
                    weights(event)*flux(event, component)
            end do
        end do
    end subroutine independent_value_oracle

end program test_fci_terminal_boundary_ledger
