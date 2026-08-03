program test_compatible_flux_elimination
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_compatible_flux_elimination, &
        assemble_compatible_flux_elimination_jvp, &
        assemble_compatible_flux_elimination_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: mass(2, 2) = reshape([ &
        2.0_dp, 0.2_dp, 0.2_dp, 1.6_dp], [2, 2])
    real(dp), parameter :: flux_to_state(2, 2) = reshape([ &
        0.8_dp, -0.2_dp, 0.1_dp, 0.6_dp], [2, 2])
    real(dp), parameter :: state_to_flux(2, 2) = reshape([ &
        0.4_dp, 0.3_dp, -0.1_dp, 0.7_dp], [2, 2])
    real(dp), parameter :: state_matrix(2, 2) = reshape([ &
        1.4_dp, 0.2_dp, 0.2_dp, 1.1_dp], [2, 2])
    real(dp), parameter :: flux_rhs(2) = [0.9_dp, -0.4_dp]
    real(dp), parameter :: state_rhs(2) = [0.2_dp, 0.7_dp]
    real(dp), parameter :: mass_dot(2, 2) = reshape([ &
        0.03_dp, -0.01_dp, -0.01_dp, 0.02_dp], [2, 2])
    real(dp), parameter :: flux_to_state_dot(2, 2) = reshape([ &
        -0.02_dp, 0.01_dp, 0.03_dp, -0.01_dp], [2, 2])
    real(dp), parameter :: state_to_flux_dot(2, 2) = reshape([ &
        0.01_dp, -0.03_dp, 0.02_dp, 0.04_dp], [2, 2])
    real(dp), parameter :: state_matrix_dot(2, 2) = reshape([ &
        -0.01_dp, 0.02_dp, 0.02_dp, -0.03_dp], [2, 2])
    real(dp), parameter :: flux_rhs_dot(2) = [0.04_dp, -0.02_dp]
    real(dp), parameter :: state_rhs_dot(2) = [-0.01_dp, 0.03_dp]
    real(dp), parameter :: recovery_bar(2) = [0.2_dp, -0.3_dp]
    real(dp), parameter :: recovery_matrix_bar(2, 2) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp], [2, 2])
    real(dp), parameter :: condensed_matrix_bar(2, 2) = reshape([ &
        -0.3_dp, 0.2_dp, 0.1_dp, 0.5_dp], [2, 2])
    real(dp), parameter :: condensed_rhs_bar(2) = [0.6_dp, -0.2_dp]
    real(dp), parameter :: eps = 1.0e-7_dp
    real(dp) :: recovery(2), recovery_matrix(2, 2), condensed_matrix(2, 2)
    real(dp) :: condensed_rhs(2), recovery_dot(2), recovery_matrix_dot(2, 2)
    real(dp) :: condensed_matrix_dot(2, 2), condensed_rhs_dot(2)
    real(dp) :: recovery_plus(2), recovery_minus(2)
    real(dp) :: recovery_matrix_plus(2, 2), recovery_matrix_minus(2, 2)
    real(dp) :: condensed_matrix_plus(2, 2), condensed_matrix_minus(2, 2)
    real(dp) :: condensed_rhs_plus(2), condensed_rhs_minus(2)
    real(dp) :: mass_bar(2, 2), flux_to_state_bar(2, 2)
    real(dp) :: state_to_flux_bar(2, 2), state_matrix_bar(2, 2)
    real(dp) :: flux_rhs_bar(2), state_rhs_bar(2)
    real(dp) :: lhs, rhs
    real(dp) :: bad_recovery(1)
    type(fortsparse_status_t) :: status

    call assemble_compatible_flux_elimination( &
        mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, state_rhs, &
        recovery, recovery_matrix, condensed_matrix, condensed_rhs, status)
    call check_condition(status%code == 0, &
        "compatible flux elimination accepts mixed blocks")
    call check_condition(maxval(abs( &
        matmul(mass, recovery) - flux_rhs)) < 1.0e-13_dp .and. &
        maxval(abs(matmul(mass, recovery_matrix) + flux_to_state)) < 1.0e-13_dp .and. &
        maxval(abs(condensed_matrix - (state_matrix + &
        matmul(state_to_flux, recovery_matrix)))) < 1.0e-13_dp .and. &
        maxval(abs(condensed_rhs - (state_rhs - &
        matmul(state_to_flux, recovery)))) < 1.0e-13_dp, &
        "compatible flux elimination matches recovery and Schur oracles")

    call assemble_compatible_flux_elimination_jvp( &
        mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, state_rhs, &
        mass_dot, flux_to_state_dot, state_to_flux_dot, state_matrix_dot, &
        flux_rhs_dot, state_rhs_dot, recovery_dot, recovery_matrix_dot, &
        condensed_matrix_dot, condensed_rhs_dot, status)
    call assemble_compatible_flux_elimination( &
        mass + eps*mass_dot, flux_to_state + eps*flux_to_state_dot, &
        state_to_flux + eps*state_to_flux_dot, state_matrix + eps*state_matrix_dot, &
        flux_rhs + eps*flux_rhs_dot, state_rhs + eps*state_rhs_dot, &
        recovery_plus, recovery_matrix_plus, condensed_matrix_plus, &
        condensed_rhs_plus, status)
    call assemble_compatible_flux_elimination( &
        mass - eps*mass_dot, flux_to_state - eps*flux_to_state_dot, &
        state_to_flux - eps*state_to_flux_dot, state_matrix - eps*state_matrix_dot, &
        flux_rhs - eps*flux_rhs_dot, state_rhs - eps*state_rhs_dot, &
        recovery_minus, recovery_matrix_minus, condensed_matrix_minus, &
        condensed_rhs_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(recovery_dot - (recovery_plus - recovery_minus)/(2.0_dp*eps))) &
        < 1.0e-8_dp .and. &
        maxval(abs(recovery_matrix_dot - (recovery_matrix_plus - &
        recovery_matrix_minus)/(2.0_dp*eps))) < 1.0e-8_dp .and. &
        maxval(abs(condensed_matrix_dot - (condensed_matrix_plus - &
        condensed_matrix_minus)/(2.0_dp*eps))) < 1.0e-8_dp .and. &
        maxval(abs(condensed_rhs_dot - (condensed_rhs_plus - &
        condensed_rhs_minus)/(2.0_dp*eps))) < 1.0e-8_dp, &
        "compatible flux elimination JVP matches central differences")

    call assemble_compatible_flux_elimination_vjp( &
        mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, state_rhs, &
        recovery, recovery_matrix, condensed_matrix, condensed_rhs, recovery_bar, &
        recovery_matrix_bar, condensed_matrix_bar, condensed_rhs_bar, mass_bar, &
        flux_to_state_bar, state_to_flux_bar, state_matrix_bar, flux_rhs_bar, &
        state_rhs_bar, status)
    lhs = sum(recovery_bar*recovery_dot) + &
        sum(recovery_matrix_bar*recovery_matrix_dot) + &
        sum(condensed_matrix_bar*condensed_matrix_dot) + &
        dot_product(condensed_rhs_bar, condensed_rhs_dot)
    rhs = sum(mass_bar*mass_dot) + sum(flux_to_state_bar*flux_to_state_dot) + &
        sum(state_to_flux_bar*state_to_flux_dot) + &
        sum(state_matrix_bar*state_matrix_dot) + &
        dot_product(flux_rhs_bar, flux_rhs_dot) + &
        dot_product(state_rhs_bar, state_rhs_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-11_dp, &
        "compatible flux elimination VJP satisfies the real dot identity")

    call assemble_compatible_flux_elimination( &
        mass, flux_to_state, state_to_flux, state_matrix, flux_rhs, state_rhs, &
        bad_recovery, recovery_matrix, condensed_matrix, condensed_rhs, status)
    call check_condition(status%code /= 0, &
        "compatible flux elimination rejects an incompatible recovery shape")
    call check_summary("Compatible flux elimination")
end program test_compatible_flux_elimination
