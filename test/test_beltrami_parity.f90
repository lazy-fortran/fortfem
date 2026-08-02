program test_beltrami_parity
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        beltrami_parity_t, evaluate_beltrami_two_region_parity, &
        validate_beltrami_parity, validate_beltrami_resonance
    use fortsparse, only: FORTSPARSE_SINGULAR, fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: nregion = 2, nsample = 3, ncomponent = 3
    integer, parameter :: nconstraint = 2
    real(dp) :: curl_hcurl(nregion, nsample, ncomponent)
    real(dp) :: curl_oracle(nregion, nsample, ncomponent)
    real(dp) :: magnetic_field(nregion, nsample, ncomponent), lambda(nregion)
    real(dp) :: divergence(nregion, nsample), divergence_target(nregion, nsample)
    real(dp) :: flux(nconstraint), flux_target(nconstraint)
    real(dp) :: helicity(nconstraint), helicity_target(nconstraint)
    real(dp) :: residual(nregion, nsample, ncomponent)
    real(dp) :: oracle_residual(nregion, nsample, ncomponent)
    real(dp) :: divergence_residual(nregion, nsample)
    real(dp) :: flux_residual(nconstraint), helicity_residual(nconstraint)
    real(dp) :: resonance_values(nregion, 2)
    type(beltrami_parity_t) :: report, invalid_report
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    call build_manufactured_state( &
        curl_hcurl, curl_oracle, magnetic_field, lambda, divergence, &
        divergence_target, flux, flux_target, helicity, helicity_target)
    call evaluate_beltrami_two_region_parity( &
        curl_hcurl, curl_oracle, magnetic_field, lambda, divergence, &
        divergence_target, flux, flux_target, helicity, helicity_target, &
        1.0e-12_dp, residual, divergence_residual, flux_residual, &
        helicity_residual, oracle_residual, report, status)
    call record_condition(status%code == 0, "two-region Beltrami parity evaluates")
    call record_condition(validate_beltrami_parity(report, status), &
        "two-region Beltrami parity report validates")
    call record_condition(report%region_count == nregion .and. &
        report%sample_count == nsample .and. report%component_count == ncomponent, &
        "parity report retains the two-region H(curl) layout")
    call record_condition(report%within_tolerance .and. &
        report%absolute_error < 1.0e-14_dp, &
        "compatible H(curl) and independent curl-eigen oracle agree")
    call record_condition(maxval(abs(residual - &
        (curl_hcurl - spread(spread(lambda, 2, nsample), 3, ncomponent) * &
        magnetic_field))) < 1.0e-14_dp, &
        "assembled residual has the independent curl-lambda-B form")
    call record_condition(maxval(abs(oracle_residual - &
        (curl_oracle - spread(spread(lambda, 2, nsample), 3, ncomponent) * &
        magnetic_field))) < 1.0e-14_dp, &
        "oracle residual is computed without reusing the assembled residual")
    call record_condition(maxval(abs(divergence_residual - &
        (divergence - divergence_target))) < 1.0e-14_dp .and. &
        maxval(abs(flux_residual - (flux - flux_target))) < 1.0e-14_dp .and. &
        maxval(abs(helicity_residual - (helicity - helicity_target))) < 1.0e-14_dp, &
        "flux, helicity, and divergence constraint rows close independently")

    curl_oracle(1, 1, 1) = curl_oracle(1, 1, 1) + 0.125_dp
    call evaluate_beltrami_two_region_parity( &
        curl_hcurl, curl_oracle, magnetic_field, lambda, divergence, &
        divergence_target, flux, flux_target, helicity, helicity_target, &
        1.0e-12_dp, residual, divergence_residual, flux_residual, &
        helicity_residual, oracle_residual, invalid_report, status)
    call record_condition(status%code == 0 .and. &
        .not. invalid_report%within_tolerance .and. &
        invalid_report%absolute_error > 0.1_dp, &
        "parity report exposes a mismatched curl path")

    resonance_values = reshape([0.25_dp, 0.8_dp, -0.7_dp, -0.2_dp], &
        shape(resonance_values))
    call validate_beltrami_resonance(lambda, resonance_values, 1.0e-10_dp, status)
    call record_condition(status%code == 0, "non-resonant region parameters are accepted")
    resonance_values(2, 1) = lambda(2)
    call validate_beltrami_resonance(lambda, resonance_values, 1.0e-10_dp, status)
    call record_condition(status%code == FORTSPARSE_SINGULAR, &
        "a supplied resonant gauge/eigenvalue is rejected")

    call check_summary("Beltrami parity")
    if (.not. all_passed) error stop 1

contains

    subroutine build_manufactured_state( &
            curl_hcurl, curl_oracle, magnetic_field, lambda, divergence, &
            divergence_target, flux, flux_target, helicity, helicity_target)
        real(dp), intent(out) :: curl_hcurl(:, :, :), curl_oracle(:, :, :)
        real(dp), intent(out) :: magnetic_field(:, :, :), lambda(:)
        real(dp), intent(out) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(out) :: flux(:), flux_target(:), helicity(:), helicity_target(:)

        curl_hcurl = reshape([ &
            0.10_dp, 0.20_dp, -0.15_dp, 0.30_dp, -0.05_dp, 0.40_dp, &
            0.25_dp, -0.35_dp, 0.45_dp, -0.20_dp, 0.15_dp, 0.05_dp, &
            -0.30_dp, 0.55_dp, 0.10_dp, 0.35_dp, -0.25_dp, 0.20_dp], &
            shape(curl_hcurl))
        magnetic_field = reshape([ &
            0.60_dp, -0.20_dp, 0.30_dp, 0.15_dp, 0.45_dp, -0.35_dp, &
            -0.25_dp, 0.50_dp, 0.40_dp, 0.30_dp, -0.10_dp, 0.20_dp, &
            0.20_dp, 0.35_dp, -0.45_dp, -0.15_dp, 0.25_dp, 0.55_dp], &
            shape(magnetic_field))
        lambda = [0.55_dp, -0.40_dp]
        curl_oracle = curl_hcurl
        divergence = reshape([0.04_dp, -0.02_dp, 0.03_dp, 0.01_dp, -0.05_dp, 0.02_dp], &
            shape(divergence))
        divergence_target = reshape([0.01_dp, -0.01_dp, 0.02_dp, 0.00_dp, -0.03_dp, 0.04_dp], &
            shape(divergence_target))
        flux = [0.9_dp, -0.4_dp]
        flux_target = [0.8_dp, -0.6_dp]
        helicity = [0.3_dp, 0.7_dp]
        helicity_target = [0.25_dp, 0.75_dp]
    end subroutine build_manufactured_state

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_beltrami_parity
