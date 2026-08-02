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

    type, public :: beltrami_shell_parity_t
        !! Geometry-labelled two-region parity and conservation ledger.
        !!
        !! The geometry label identifies the caller's sample set only.  This
        !! contract does not construct slab or toroidal-shell coordinates.
        character(len=32) :: schema_version = "fortfem-beltrami-shell-1"
        character(len=32) :: geometry_kind = ""
        integer :: region_count = 0
        integer :: sample_count = 0
        integer :: component_count = 0
        real(dp) :: tolerance = 0.0_dp
        real(dp) :: weighted_hcurl_residual_norm = 0.0_dp
        real(dp) :: weighted_oracle_residual_norm = 0.0_dp
        real(dp) :: weighted_absolute_error = 0.0_dp
        real(dp) :: weighted_relative_error = 0.0_dp
        real(dp) :: flux_closure_norm = 0.0_dp
        real(dp) :: helicity_closure_norm = 0.0_dp
        real(dp) :: energy_closure_norm = 0.0_dp
        real(dp) :: ledger_closure_norm = 0.0_dp
        real(dp) :: hcurl_energy = 0.0_dp
        real(dp) :: oracle_energy = 0.0_dp
        real(dp) :: energy_parity_error = 0.0_dp
        real(dp), allocatable :: energy_by_region(:)
        logical :: ledger_closed = .false.
        logical :: within_tolerance = .false.
    end type beltrami_shell_parity_t

    public :: compare_beltrami_two_region_residual
    public :: validate_beltrami_parity
    public :: validate_beltrami_resonance
    public :: compare_beltrami_shell_residual
    public :: validate_beltrami_shell_parity

contains

    subroutine compare_beltrami_two_region_residual( &
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
    end subroutine compare_beltrami_two_region_residual

    subroutine compare_beltrami_shell_residual( &
            geometry_kind, curl_hcurl, curl_oracle, magnetic_field, lambda, &
            sample_weight, divergence, divergence_target, flux, flux_target, &
            helicity, helicity_target, energy_target, tolerance, report, status)
        !! Compare compatible and independent curl paths on caller-supplied
        !! slab or toroidal-shell samples and close a neutral ledger.
        !!
        !! `energy_by_region` is the weighted quadratic field energy
        !! 1/2 sum_q w_q B_q.B_q.  Flux, helicity, and energy targets are
        !! supplied by the consuming client; no coordinate or equilibrium
        !! model is inferred here.
        character(len=*), intent(in) :: geometry_kind
        real(dp), intent(in) :: curl_hcurl(:, :, :), curl_oracle(:, :, :)
        real(dp), intent(in) :: magnetic_field(:, :, :), lambda(:)
        real(dp), intent(in) :: sample_weight(:, :)
        real(dp), intent(in) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(in) :: flux(:), flux_target(:), helicity(:), helicity_target(:)
        real(dp), intent(in) :: energy_target(:), tolerance
        type(beltrami_shell_parity_t), intent(out) :: report
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: residual(:, :, :), oracle_residual(:, :, :)
        real(dp), allocatable :: divergence_residual(:, :), flux_residual(:)
        real(dp), allocatable :: helicity_residual(:)
        real(dp) :: hcurl_norm_sq, oracle_norm_sq, difference_norm_sq
        real(dp) :: energy_error_sq
        type(beltrami_parity_t) :: base_report
        integer :: region, sample, component

        report = beltrami_shell_parity_t()
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Beltrami shell parity received incompatible arrays")
        if (trim(geometry_kind) /= "slab" .and. &
                trim(geometry_kind) /= "toroidal-shell") return
        if (.not. valid_shell_shapes(curl_hcurl, curl_oracle, magnetic_field, lambda, &
                sample_weight, divergence, divergence_target, flux, flux_target, &
                helicity, helicity_target, energy_target)) return
        if (.not. finite_shell_inputs(curl_hcurl, curl_oracle, magnetic_field, lambda, &
                sample_weight, divergence, divergence_target, flux, flux_target, &
                helicity, helicity_target, energy_target) .or. &
                .not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp .or. &
                any(sample_weight <= 0.0_dp)) return

        allocate(residual(size(curl_hcurl, 1), size(curl_hcurl, 2), &
            size(curl_hcurl, 3)))
        allocate(oracle_residual(size(curl_hcurl, 1), size(curl_hcurl, 2), &
            size(curl_hcurl, 3)))
        allocate(divergence_residual(size(divergence, 1), size(divergence, 2)))
        allocate(flux_residual(size(flux)), helicity_residual(size(helicity)))
        call compare_beltrami_two_region_residual( &
            curl_hcurl, curl_oracle, magnetic_field, lambda, divergence, &
            divergence_target, flux, flux_target, helicity, helicity_target, &
            tolerance, &
            residual, divergence_residual, flux_residual, helicity_residual, &
            oracle_residual, base_report, status)
        if (status%code /= FORTSPARSE_OK) return

        hcurl_norm_sq = 0.0_dp
        oracle_norm_sq = 0.0_dp
        difference_norm_sq = 0.0_dp
        if (.not. allocated(report%energy_by_region)) allocate( &
            report%energy_by_region(size(magnetic_field, 1)))
        report%energy_by_region = 0.0_dp
        do region = 1, size(curl_hcurl, 1)
            do sample = 1, size(curl_hcurl, 2)
                do component = 1, size(curl_hcurl, 3)
                    hcurl_norm_sq = hcurl_norm_sq + sample_weight(region, sample) * &
                        residual(region, sample, component)**2
                    oracle_norm_sq = oracle_norm_sq + sample_weight(region, sample) * &
                        oracle_residual(region, sample, component)**2
                    difference_norm_sq = difference_norm_sq + &
                        sample_weight(region, sample) * &
                        (residual(region, sample, component) - &
                        oracle_residual(region, sample, component))**2
                    report%energy_by_region(region) = &
                        report%energy_by_region(region) + &
                        0.5_dp * sample_weight(region, sample) * &
                        magnetic_field(region, sample, component)**2
                end do
            end do
        end do
        energy_error_sq = sum((report%energy_by_region - energy_target)**2)
        report%schema_version = "fortfem-beltrami-shell-1"
        report%geometry_kind = trim(geometry_kind)
        report%region_count = size(curl_hcurl, 1)
        report%sample_count = size(curl_hcurl, 2)
        report%component_count = size(curl_hcurl, 3)
        report%tolerance = tolerance
        report%weighted_hcurl_residual_norm = sqrt(hcurl_norm_sq)
        report%weighted_oracle_residual_norm = sqrt(oracle_norm_sq)
        report%weighted_absolute_error = sqrt(difference_norm_sq)
        report%weighted_relative_error = report%weighted_absolute_error / &
            max(report%weighted_oracle_residual_norm, tiny(1.0_dp))
        report%flux_closure_norm = sqrt(sum((flux - flux_target)**2))
        report%helicity_closure_norm = sqrt(sum((helicity - helicity_target)**2))
        report%energy_closure_norm = sqrt(energy_error_sq)
        report%ledger_closure_norm = sqrt(report%flux_closure_norm**2 + &
            report%helicity_closure_norm**2 + report%energy_closure_norm**2)
        report%hcurl_energy = sum(report%energy_by_region)
        report%oracle_energy = report%hcurl_energy
        report%energy_parity_error = abs(report%hcurl_energy - report%oracle_energy)
        report%ledger_closed = report%ledger_closure_norm <= tolerance
        report%within_tolerance = report%ledger_closed .and. &
            (report%weighted_absolute_error <= tolerance .or. &
            report%weighted_relative_error <= tolerance)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compare_beltrami_shell_residual

    logical function validate_beltrami_shell_parity(report, status) result(valid)
        type(beltrami_shell_parity_t), intent(in) :: report
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: values(12)

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Beltrami shell parity report is invalid")
        if (report%schema_version /= "fortfem-beltrami-shell-1" .or. &
                (trim(report%geometry_kind) /= "slab" .and. &
                trim(report%geometry_kind) /= "toroidal-shell") .or. &
                report%region_count /= 2 .or. report%sample_count < 1 .or. &
                report%component_count < 1 .or. report%tolerance < 0.0_dp) return
        if (.not. allocated(report%energy_by_region)) return
        if (size(report%energy_by_region) /= report%region_count) return
        values = [report%tolerance, report%weighted_hcurl_residual_norm, &
            report%weighted_oracle_residual_norm, report%weighted_absolute_error, &
            report%weighted_relative_error, report%flux_closure_norm, &
            report%helicity_closure_norm, report%energy_closure_norm, &
            report%ledger_closure_norm, report%hcurl_energy, report%oracle_energy, &
            report%energy_parity_error]
        if (.not. all(ieee_is_finite(values)) .or. &
                .not. all(ieee_is_finite(report%energy_by_region)) .or. &
                any(values(2:) < 0.0_dp) .or. &
                any(report%energy_by_region < 0.0_dp)) return
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_beltrami_shell_parity

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

    logical function valid_shell_shapes( &
            curl_hcurl, curl_oracle, magnetic_field, lambda, sample_weight, divergence, &
            divergence_target, flux, flux_target, helicity, helicity_target, &
            energy_target) result(valid)
        real(dp), intent(in) :: curl_hcurl(:, :, :), curl_oracle(:, :, :)
        real(dp), intent(in) :: magnetic_field(:, :, :), lambda(:)
        real(dp), intent(in) :: sample_weight(:, :)
        real(dp), intent(in) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(in) :: flux(:), flux_target(:), helicity(:), helicity_target(:)
        real(dp), intent(in) :: energy_target(:)

        valid = .false.
        if (size(curl_hcurl, 1) /= 2) return
        if (size(curl_hcurl, 2) < 1 .or. size(curl_hcurl, 3) < 1) return
        if (.not. all(shape(curl_oracle) == shape(curl_hcurl))) return
        if (.not. all(shape(magnetic_field) == shape(curl_hcurl))) return
        if (.not. all(shape(sample_weight) == shape(curl_hcurl(:, :, 1)))) return
        if (size(lambda) /= size(curl_hcurl, 1)) return
        if (.not. all(shape(divergence) == shape(divergence_target))) return
        if (size(divergence, 1) /= size(curl_hcurl, 1) .or. &
                size(divergence, 2) /= size(curl_hcurl, 2)) return
        if (size(flux) /= size(flux_target) .or. &
                size(helicity) /= size(helicity_target)) return
        if (size(energy_target) /= size(curl_hcurl, 1)) return
        valid = .true.
    end function valid_shell_shapes

    logical function finite_shell_inputs( &
            curl_hcurl, curl_oracle, magnetic_field, lambda, sample_weight, divergence, &
            divergence_target, flux, flux_target, helicity, helicity_target, &
            energy_target) result(valid)
        real(dp), intent(in) :: curl_hcurl(:, :, :), curl_oracle(:, :, :)
        real(dp), intent(in) :: magnetic_field(:, :, :), lambda(:)
        real(dp), intent(in) :: sample_weight(:, :)
        real(dp), intent(in) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(in) :: flux(:), flux_target(:), helicity(:), helicity_target(:)
        real(dp), intent(in) :: energy_target(:)

        valid = finite_inputs(curl_hcurl, curl_oracle, magnetic_field, lambda, divergence, &
            divergence_target, flux, flux_target, helicity, helicity_target) .and. &
            all(ieee_is_finite(sample_weight)) .and. all(ieee_is_finite(energy_target))
    end function finite_shell_inputs

    logical function finite_scalars(values) result(valid)
        real(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(values))
    end function finite_scalars

end module fortfem_beltrami_parity
