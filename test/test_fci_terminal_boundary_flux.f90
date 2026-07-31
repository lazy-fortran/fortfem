program test_fci_terminal_boundary_flux
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_fci_terminal_boundary_flux, &
        assemble_fci_terminal_boundary_flux_jvp, &
        assemble_fci_terminal_boundary_flux_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: owners(4) = [1, 2, 1, 3]
    real(dp), parameter :: weights(4) = [0.5_dp, 1.2_dp, 0.25_dp, 0.75_dp]
    real(dp), parameter :: terminal_flux(4) = [2.0_dp, -1.0_dp, 0.5_dp, 3.0_dp]
    real(dp), parameter :: canonical_volumes(3) = [2.0_dp, 4.0_dp, 1.5_dp]
    real(dp), parameter :: weights_dot(4) = [0.03_dp, -0.04_dp, 0.02_dp, 0.01_dp]
    real(dp), parameter :: terminal_flux_dot(4) = [-0.2_dp, 0.1_dp, 0.05_dp, -0.3_dp]
    real(dp), parameter :: canonical_volumes_dot(3) = [0.1_dp, -0.2_dp, 0.07_dp]
    real(dp), parameter :: contribution_bar(3) = [0.4_dp, -0.6_dp, 0.8_dp]
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: contribution(3), contribution_dot(3), contribution_plus(3)
    real(dp) :: contribution_minus(3), weights_bar(4), terminal_flux_bar(4)
    real(dp) :: canonical_volumes_bar(3), left, right
    real(dp) :: weighted_integral, contribution_dot_integral
    type(fortsparse_status_t) :: status

    call assemble_fci_terminal_boundary_flux( &
        owners, weights, terminal_flux, canonical_volumes, contribution, status)
    call check_condition(status%code == 0, &
        "FCI terminal flux assembly accepts fixed owner indices")
    call check_condition(maxval(abs(contribution - &
        [0.5625_dp, -0.3_dp, 1.5_dp])) < 1.0e-14_dp, &
        "FCI terminal flux assembly matches the independent volume-scaled oracle")
    weighted_integral = dot_product(canonical_volumes, contribution)
    call check_condition(abs(weighted_integral - dot_product(weights, terminal_flux)) < &
        1.0e-14_dp, "FCI terminal flux assembly preserves the weighted balance")

    call assemble_fci_terminal_boundary_flux_jvp( &
        owners, weights, terminal_flux, canonical_volumes, weights_dot, &
        terminal_flux_dot, canonical_volumes_dot, contribution_dot, status)
    call check_condition(status%code == 0, &
        "FCI terminal flux JVP accepts fixed owner topology")
    call assemble_fci_terminal_boundary_flux( &
        owners, weights + step*weights_dot, terminal_flux + step*terminal_flux_dot, &
        canonical_volumes + step*canonical_volumes_dot, contribution_plus, status)
    call assemble_fci_terminal_boundary_flux( &
        owners, weights - step*weights_dot, terminal_flux - step*terminal_flux_dot, &
        canonical_volumes - step*canonical_volumes_dot, contribution_minus, status)
    call check_condition(maxval(abs(contribution_dot - (contribution_plus - &
        contribution_minus)/(2.0_dp*step))) < 5.0e-9_dp, &
        "FCI terminal flux JVP matches an independent central difference")
    contribution_dot_integral = dot_product(canonical_volumes, contribution_dot) + &
        dot_product(canonical_volumes_dot, contribution)
    call check_condition(abs(contribution_dot_integral - &
        (dot_product(weights_dot, terminal_flux) + &
        dot_product(weights, terminal_flux_dot))) < 1.0e-12_dp, &
        "FCI terminal flux JVP preserves the differentiated weighted balance")

    call assemble_fci_terminal_boundary_flux_vjp( &
        owners, weights, terminal_flux, canonical_volumes, contribution_bar, &
        weights_bar, terminal_flux_bar, canonical_volumes_bar, status)
    left = dot_product(contribution_bar, contribution_dot)
    right = dot_product(weights_bar, weights_dot) + &
        dot_product(terminal_flux_bar, terminal_flux_dot) + &
        dot_product(canonical_volumes_bar, canonical_volumes_dot)
    call check_condition(status%code == 0 .and. abs(left - right) < 1.0e-12_dp, &
        "FCI terminal flux VJP satisfies the real dot-product identity")

    call assemble_fci_terminal_boundary_flux( &
        [1, 4, 1, 2], weights, terminal_flux, canonical_volumes, contribution, status)
    call check_condition(status%code /= 0, &
        "FCI terminal flux assembly rejects an invalid owner")
    call assemble_fci_terminal_boundary_flux( &
        owners, [-0.1_dp, 1.2_dp, 0.25_dp, 0.75_dp], terminal_flux, &
        canonical_volumes, contribution, status)
    call check_condition(status%code /= 0, &
        "FCI terminal flux assembly rejects a non-positive weight")
    call check_summary("FCI terminal boundary flux")
end program test_fci_terminal_boundary_flux
