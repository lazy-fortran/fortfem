program test_beltrami_facade
    !! Smoke-test the canonical FEEC facade for the Beltrami parity contracts.
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        beltrami_parity_t, beltrami_shell_parity_t, &
        compare_beltrami_two_region_residual, &
        compare_beltrami_shell_residual, &
        validate_beltrami_parity, validate_beltrami_shell_parity
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp) :: curl_hcurl(2, 1, 3), curl_oracle(2, 1, 3)
    real(dp) :: magnetic_field(2, 1, 3), lambda(2)
    real(dp) :: divergence(2, 1), divergence_target(2, 1)
    real(dp) :: flux(1), flux_target(1), helicity(1), helicity_target(1)
    real(dp) :: sample_weight(2, 1), energy_target(2)
    real(dp) :: residual(2, 1, 3), oracle_residual(2, 1, 3)
    real(dp) :: divergence_residual(2, 1), flux_residual(1), helicity_residual(1)
    type(beltrami_parity_t) :: parity_report
    type(beltrami_shell_parity_t) :: shell_report
    type(fortsparse_status_t) :: status

    magnetic_field(1, 1, :) = [1.0_dp, -0.5_dp, 0.25_dp]
    magnetic_field(2, 1, :) = [0.5_dp, 0.25_dp, -1.0_dp]
    lambda = [2.0_dp, -1.0_dp]
    curl_hcurl(1, 1, :) = lambda(1) * magnetic_field(1, 1, :)
    curl_hcurl(2, 1, :) = lambda(2) * magnetic_field(2, 1, :)
    curl_oracle = curl_hcurl
    divergence = 0.0_dp
    divergence_target = divergence
    flux = [0.75_dp]
    flux_target = flux
    helicity = [-0.25_dp]
    helicity_target = helicity

    call compare_beltrami_two_region_residual( &
        curl_hcurl, curl_oracle, magnetic_field, lambda, divergence, &
        divergence_target, flux, flux_target, helicity, helicity_target, &
        1.0e-12_dp, residual, divergence_residual, flux_residual, &
        helicity_residual, oracle_residual, parity_report, status)
    call check_condition(status%code == 0 .and. &
        validate_beltrami_parity(parity_report, status) .and. &
        parity_report%within_tolerance, &
        "canonical FEEC facade exports the two-region Beltrami comparator")

    sample_weight = 1.0_dp
    energy_target = [0.5_dp * sum(magnetic_field(1, 1, :)**2), &
        0.5_dp * sum(magnetic_field(2, 1, :)**2)]
    call compare_beltrami_shell_residual( &
        "slab", curl_hcurl, curl_oracle, magnetic_field, lambda, sample_weight, &
        divergence, divergence_target, flux, flux_target, helicity, &
        helicity_target, energy_target, 1.0e-12_dp, shell_report, status)
    call check_condition(status%code == 0 .and. &
        validate_beltrami_shell_parity(shell_report, status) .and. &
        shell_report%ledger_closed, &
        "canonical FEEC facade exports the shell energy ledger")

    call check_summary("Beltrami FEEC facade")
end program test_beltrami_facade
