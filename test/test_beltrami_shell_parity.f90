program test_beltrami_shell_parity
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        beltrami_shell_parity_t, evaluate_beltrami_shell_parity, &
        validate_beltrami_shell_parity
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: nregion = 2, nsample = 4, ncomponent = 3
    integer, parameter :: nconstraint = 2
    real(dp) :: curl_hcurl(nregion, nsample, ncomponent)
    real(dp) :: curl_oracle(nregion, nsample, ncomponent)
    real(dp) :: magnetic_field(nregion, nsample, ncomponent), lambda(nregion)
    real(dp) :: sample_weight(nregion, nsample)
    real(dp) :: divergence(nregion, nsample), divergence_target(nregion, nsample)
    real(dp) :: flux(nconstraint), flux_target(nconstraint)
    real(dp) :: helicity(nconstraint), helicity_target(nconstraint)
    real(dp) :: energy_target(nregion)
    real(dp) :: expected_energy(nregion)
    type(beltrami_shell_parity_t) :: report
    type(fortsparse_status_t) :: status

    call build_state(curl_hcurl, curl_oracle, magnetic_field, lambda, sample_weight, &
        divergence, divergence_target, flux, flux_target, helicity, helicity_target, &
        energy_target, expected_energy)

    call evaluate_beltrami_shell_parity( &
        "slab", curl_hcurl, curl_oracle, magnetic_field, lambda, sample_weight, &
        divergence, divergence_target, flux, flux_target, helicity, helicity_target, &
        energy_target, 1.0e-12_dp, report, status)
    call check_condition(status%code == 0, "two-region slab ledger evaluates")
    call check_condition(validate_beltrami_shell_parity(report, status), &
        "two-region slab ledger validates")
    call check_condition(report%geometry_kind == "slab" .and. &
        report%region_count == 2, &
        "slab report retains the two-region geometry contract")
    call check_condition(report%within_tolerance .and. report%ledger_closed, &
        "compatible H(curl) and curl-eigen paths close the slab ledger")
    call check_condition(maxval(abs(report%energy_by_region - expected_energy)) < &
        1.0e-14_dp, &
        "independent magnetic-energy quadrature closes each slab region")

    call evaluate_beltrami_shell_parity( &
        "toroidal-shell", curl_hcurl, curl_oracle, magnetic_field, lambda, &
        sample_weight, &
        divergence, divergence_target, flux, flux_target, helicity, helicity_target, &
        energy_target, 1.0e-12_dp, report, status)
    call check_condition(status%code == 0 .and. &
        report%geometry_kind == "toroidal-shell", &
        "toroidal-shell ledger uses the same neutral parity contract")
    call check_condition(report%within_tolerance .and. &
        report%energy_parity_error < 1.0e-14_dp, &
        "toroidal-shell compatible and independent energies agree")

    curl_oracle(2, 3, 2) = curl_oracle(2, 3, 2) + 0.125_dp
    call evaluate_beltrami_shell_parity( &
        "toroidal-shell", curl_hcurl, curl_oracle, magnetic_field, lambda, &
        sample_weight, &
        divergence, divergence_target, flux, flux_target, helicity, helicity_target, &
        energy_target, 1.0e-12_dp, report, status)
    call check_condition(status%code == 0 .and. .not. report%within_tolerance .and. &
        report%weighted_absolute_error > 0.09_dp, &
        "toroidal-shell report exposes an independent curl mismatch")

    energy_target(1) = energy_target(1) + 0.25_dp
    call evaluate_beltrami_shell_parity( &
        "slab", curl_hcurl, curl_hcurl, magnetic_field, lambda, sample_weight, &
        divergence, divergence_target, flux, flux_target, helicity, helicity_target, &
        energy_target, 1.0e-12_dp, report, status)
    call check_condition(status%code == 0 .and. .not. report%ledger_closed .and. &
        report%energy_closure_norm > 0.2_dp, &
        "slab report rejects an energy-closure mismatch independently")

    call check_summary("Beltrami shell parity")

contains

    subroutine build_state( &
            curl_hcurl, curl_oracle, magnetic_field, lambda, sample_weight, &
            divergence, &
            divergence_target, flux, flux_target, helicity, helicity_target, &
            energy_target, expected_energy)
        real(dp), intent(out) :: curl_hcurl(:, :, :), curl_oracle(:, :, :)
        real(dp), intent(out) :: magnetic_field(:, :, :), lambda(:)
        real(dp), intent(out) :: sample_weight(:, :)
        real(dp), intent(out) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(out) :: flux(:), flux_target(:), helicity(:), &
            helicity_target(:)
        real(dp), intent(out) :: energy_target(:), expected_energy(:)
        integer :: region, sample

        magnetic_field = reshape([ &
            1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
            1.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, &
            0.5_dp, -0.5_dp, 1.0_dp, 0.0_dp, 2.0_dp, -1.0_dp, &
            1.0_dp, 0.0_dp, 2.0_dp, -0.5_dp, 0.5_dp, 0.5_dp], &
            shape(magnetic_field))
        lambda = [1.5_dp, -0.75_dp]
        do region = 1, nregion
            do sample = 1, nsample
                curl_hcurl(region, sample, :) = lambda(region) * &
                    magnetic_field(region, sample, :)
            end do
        end do
        curl_oracle = curl_hcurl
        sample_weight = reshape([ &
            0.25_dp, 0.50_dp, 0.75_dp, 1.00_dp, &
            0.40_dp, 0.60_dp, 0.80_dp, 1.20_dp], shape(sample_weight))
        divergence = 0.0_dp
        divergence_target = 0.0_dp
        flux = [1.25_dp, -0.5_dp]
        flux_target = flux
        helicity = [0.75_dp, 1.5_dp]
        helicity_target = helicity
        do region = 1, nregion
            expected_energy(region) = 0.0_dp
            do sample = 1, nsample
                expected_energy(region) = expected_energy(region) + 0.5_dp * &
                    sample_weight(region, sample) * &
                    sum(magnetic_field(region, sample, :)**2)
            end do
        end do
        energy_target = expected_energy
    end subroutine build_state

end program test_beltrami_shell_parity
