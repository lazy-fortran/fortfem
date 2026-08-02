module fortfem_beltrami_parity
    !! Two-path manufactured parity diagnostics for multi-region Beltrami blocks.
    !!
    !! The compatible H(curl) path supplies `curl_hcurl` to the neutral
    !! `assemble_beltrami_residual` contract.  A second, independent algebraic
    !! path evaluates `curl_oracle-lambda*B` directly.  This module compares
    !! those paths and reports the supplied constraint-row norms; it does not
    !! construct a field, gauge, topology, or physical equilibrium.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_beltrami_residual, only: assemble_beltrami_residual
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        FORTSPARSE_SINGULAR, fortsparse_status_t, status_set
    implicit none
    private

    type, public :: beltrami_parity_t
        character(len=32) :: schema_version = "fortfem-beltrami-parity-1"
        integer :: region_count = 0
        integer :: sample_count = 0
        integer :: component_count = 0
        real(dp) :: tolerance = 0.0_dp
        real(dp) :: hcurl_residual_norm = 0.0_dp
        real(dp) :: oracle_residual_norm = 0.0_dp
        real(dp) :: absolute_error = 0.0_dp
        real(dp) :: relative_error = 0.0_dp
        real(dp) :: divergence_residual_norm = 0.0_dp
        real(dp) :: flux_residual_norm = 0.0_dp
        real(dp) :: helicity_residual_norm = 0.0_dp
        logical :: within_tolerance = .false.
    end type beltrami_parity_t

    public :: evaluate_beltrami_two_region_parity
    public :: validate_beltrami_parity
    public :: validate_beltrami_resonance

contains

    subroutine evaluate_beltrami_two_region_parity( &
            curl_hcurl, curl_oracle, magnetic_field, lambda, divergence, &
            divergence_target, flux, flux_target, helicity, helicity_target, &
            tolerance, residual, divergence_residual, flux_residual, &
            helicity_residual, oracle_residual, report, status)
        real(dp), intent(in) :: curl_hcurl(:, :, :), curl_oracle(:, :, :)
        real(dp), intent(in) :: magnetic_field(:, :, :), lambda(:)
        real(dp), intent(in) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(in) :: flux(:), flux_target(:), helicity(:), helicity_target(:)
        real(dp), intent(in) :: tolerance
        real(dp), intent(out) :: residual(:, :, :), divergence_residual(:, :)
        real(dp), intent(out) :: flux_residual(:), helicity_residual(:)
        real(dp), intent(out) :: oracle_residual(:, :, :)
        type(beltrami_parity_t), intent(out) :: report
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: difference_norm
        integer :: region, sample, component

        report = beltrami_parity_t()
        residual = 0.0_dp
        divergence_residual = 0.0_dp
        flux_residual = 0.0_dp
        helicity_residual = 0.0_dp
        oracle_residual = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Beltrami parity received incompatible arrays")
        if (.not. valid_shapes(curl_hcurl, curl_oracle, magnetic_field, lambda, &
                divergence, divergence_target, flux, flux_target, helicity, &
                helicity_target, residual, divergence_residual, flux_residual, &
                helicity_residual, oracle_residual)) return
        if (.not. finite_inputs(curl_hcurl, curl_oracle, magnetic_field, lambda, &
                divergence, divergence_target, flux, flux_target, helicity, &
                helicity_target) .or. .not. ieee_is_finite(tolerance) .or. &
                tolerance < 0.0_dp) return

        call assemble_beltrami_residual( &
            curl_hcurl, magnetic_field, lambda, divergence, divergence_target, &
            flux, flux_target, helicity, helicity_target, residual, &
            divergence_residual, flux_residual, helicity_residual, status)
        if (status%code /= FORTSPARSE_OK) return

        do region = 1, size(curl_hcurl, 1)
            do sample = 1, size(curl_hcurl, 2)
                do component = 1, size(curl_hcurl, 3)
                    oracle_residual(region, sample, component) = &
                        curl_oracle(region, sample, component) - &
                        lambda(region)*magnetic_field(region, sample, component)
                end do
            end do
        end do
        difference_norm = sqrt(sum((residual - oracle_residual)**2))
        report%region_count = size(curl_hcurl, 1)
        report%sample_count = size(curl_hcurl, 2)
        report%component_count = size(curl_hcurl, 3)
        report%tolerance = tolerance
        report%hcurl_residual_norm = sqrt(sum(residual**2))
        report%oracle_residual_norm = sqrt(sum(oracle_residual**2))
        report%absolute_error = difference_norm
        report%relative_error = difference_norm/max(report%oracle_residual_norm, tiny(1.0_dp))
        report%divergence_residual_norm = sqrt(sum(divergence_residual**2))
        report%flux_residual_norm = sqrt(sum(flux_residual**2))
        report%helicity_residual_norm = sqrt(sum(helicity_residual**2))
        report%within_tolerance = report%absolute_error <= tolerance .or. &
            report%relative_error <= tolerance
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_beltrami_two_region_parity

    logical function validate_beltrami_parity(report, status) result(valid)
        type(beltrami_parity_t), intent(in) :: report
        type(fortsparse_status_t), intent(out) :: status

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Beltrami parity report is invalid")
        if (report%schema_version /= "fortfem-beltrami-parity-1") return
        if (report%region_count < 1 .or. report%sample_count < 1 .or. &
                report%component_count < 1 .or. report%tolerance < 0.0_dp) return
        if (.not. finite_scalars([report%tolerance, report%hcurl_residual_norm, &
                report%oracle_residual_norm, report%absolute_error, report%relative_error, &
                report%divergence_residual_norm, report%flux_residual_norm, &
                report%helicity_residual_norm])) return
        if (report%hcurl_residual_norm < 0.0_dp .or. report%oracle_residual_norm < 0.0_dp .or. &
                report%absolute_error < 0.0_dp .or. report%relative_error < 0.0_dp .or. &
                report%divergence_residual_norm < 0.0_dp .or. &
                report%flux_residual_norm < 0.0_dp .or. report%helicity_residual_norm < 0.0_dp) return
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_beltrami_parity

    subroutine validate_beltrami_resonance(lambda, resonance_values, tolerance, status)
        !! Reject a supplied region parameter within `tolerance` of a forbidden
        !! gauge/eigenvalue.  The caller supplies the spectrum; no gauge choice
        !! or eigenproblem is inferred here.
        real(dp), intent(in) :: lambda(:), resonance_values(:, :), tolerance
        type(fortsparse_status_t), intent(out) :: status
        integer :: region, mode

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Beltrami resonance metadata is incompatible")
        if (size(lambda) < 1 .or. size(resonance_values, 1) /= size(lambda) .or. &
                size(resonance_values, 2) < 1 .or. .not. ieee_is_finite(tolerance) .or. &
                tolerance < 0.0_dp) return
        if (.not. all(ieee_is_finite(lambda)) .or. &
                .not. all(ieee_is_finite(resonance_values))) return
        do region = 1, size(lambda)
            do mode = 1, size(resonance_values, 2)
                if (abs(lambda(region) - resonance_values(region, mode)) <= tolerance) then
                    call status_set(status, FORTSPARSE_SINGULAR, &
                        "Beltrami parameter is resonant with supplied gauge/eigenvalue")
                    return
                end if
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_beltrami_resonance

    logical function valid_shapes( &
            curl_hcurl, curl_oracle, magnetic_field, lambda, divergence, &
            divergence_target, flux, flux_target, helicity, helicity_target, &
            residual, divergence_residual, flux_residual, helicity_residual, &
            oracle_residual) result(valid)
        real(dp), intent(in) :: curl_hcurl(:, :, :), curl_oracle(:, :, :)
        real(dp), intent(in) :: magnetic_field(:, :, :), lambda(:)
        real(dp), intent(in) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(in) :: flux(:), flux_target(:), helicity(:), helicity_target(:)
        real(dp), intent(in) :: residual(:, :, :), divergence_residual(:, :)
        real(dp), intent(in) :: flux_residual(:), helicity_residual(:)
        real(dp), intent(in) :: oracle_residual(:, :, :)

        valid = size(curl_hcurl, 1) > 0 .and. size(curl_hcurl, 2) > 0 .and. &
            size(curl_hcurl, 3) > 0 .and. all(shape(curl_oracle) == shape(curl_hcurl)) .and. &
            all(shape(magnetic_field) == shape(curl_hcurl)) .and. &
            all(shape(residual) == shape(curl_hcurl)) .and. &
            all(shape(oracle_residual) == shape(curl_hcurl)) .and. &
            size(lambda) == size(curl_hcurl, 1) .and. &
            all(shape(divergence_target) == shape(divergence)) .and. &
            all(shape(divergence_residual) == shape(divergence)) .and. &
            size(flux_target) == size(flux) .and. size(flux_residual) == size(flux) .and. &
            size(helicity_target) == size(helicity) .and. &
            size(helicity_residual) == size(helicity)
    end function valid_shapes

    logical function finite_inputs( &
            curl_hcurl, curl_oracle, magnetic_field, lambda, divergence, &
            divergence_target, flux, flux_target, helicity, helicity_target) result(valid)
        real(dp), intent(in) :: curl_hcurl(:, :, :), curl_oracle(:, :, :)
        real(dp), intent(in) :: magnetic_field(:, :, :), lambda(:)
        real(dp), intent(in) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(in) :: flux(:), flux_target(:), helicity(:), helicity_target(:)

        valid = all(ieee_is_finite(curl_hcurl)) .and. all(ieee_is_finite(curl_oracle)) .and. &
            all(ieee_is_finite(magnetic_field)) .and. all(ieee_is_finite(lambda)) .and. &
            all(ieee_is_finite(divergence)) .and. all(ieee_is_finite(divergence_target)) .and. &
            all(ieee_is_finite(flux)) .and. all(ieee_is_finite(flux_target)) .and. &
            all(ieee_is_finite(helicity)) .and. all(ieee_is_finite(helicity_target))
    end function finite_inputs

    logical function finite_scalars(values) result(valid)
        real(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(values))
    end function finite_scalars

end module fortfem_beltrami_parity
