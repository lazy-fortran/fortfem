module fortfem_xfem_blending_correction
    !! Fixed-topology corrected XFEM blending enrichment.
    !!
    !! A partition-of-unity ramp is the sum of base functions whose nodes are
    !! enriched.  It is zero in a standard element, one in a fully enriched
    !! element, and varies continuously in a blending element.  The discrete
    !! enrichment mask is caller-owned and is not differentiated.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_blending_corrected_enrichment
    public :: evaluate_blending_corrected_enrichment_jvp
    public :: evaluate_blending_corrected_enrichment_vjp

contains

    subroutine evaluate_blending_corrected_enrichment( &
            base_values, enriched_mask, enrichment_values, corrected_values, &
            status)
        real(dp), intent(in) :: base_values(:, :), enrichment_values(:)
        logical, intent(in) :: enriched_mask(:)
        real(dp), intent(out) :: corrected_values(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: ramp(:)

        corrected_values = 0.0_dp
        call validate_blending_inputs( &
            base_values, enriched_mask, enrichment_values, corrected_values, status)
        if (status%code /= FORTSPARSE_OK) return
        allocate(ramp(size(enrichment_values)))
        call blending_ramp(base_values, enriched_mask, ramp)
        corrected_values = ramp*enrichment_values
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_blending_corrected_enrichment

    subroutine evaluate_blending_corrected_enrichment_jvp( &
            base_values, enriched_mask, enrichment_values, base_dot, &
            enrichment_dot, corrected_dot, status)
        real(dp), intent(in) :: base_values(:, :), enrichment_values(:)
        logical, intent(in) :: enriched_mask(:)
        real(dp), intent(in) :: base_dot(:, :), enrichment_dot(:)
        real(dp), intent(out) :: corrected_dot(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: ramp(:), ramp_dot(:)

        corrected_dot = 0.0_dp
        call validate_blending_inputs( &
            base_values, enriched_mask, enrichment_values, corrected_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(base_dot, 1) /= size(base_values, 1) .or. &
            size(base_dot, 2) /= size(base_values, 2) .or. &
            size(enrichment_dot) /= size(enrichment_values) .or. &
            any(.not. ieee_is_finite(base_dot)) .or. &
            any(.not. ieee_is_finite(enrichment_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "corrected enrichment JVP has incompatible increments")
            return
        end if
        allocate(ramp(size(enrichment_values)), ramp_dot(size(enrichment_values)))
        call blending_ramp(base_values, enriched_mask, ramp)
        call blending_ramp(base_dot, enriched_mask, ramp_dot)
        corrected_dot = ramp_dot*enrichment_values + ramp*enrichment_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_blending_corrected_enrichment_jvp

    subroutine evaluate_blending_corrected_enrichment_vjp( &
            base_values, enriched_mask, enrichment_values, corrected_bar, &
            base_bar, enrichment_bar, status)
        real(dp), intent(in) :: base_values(:, :), enrichment_values(:)
        logical, intent(in) :: enriched_mask(:)
        real(dp), intent(in) :: corrected_bar(:)
        real(dp), intent(out) :: base_bar(:, :), enrichment_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: ramp(:)
        integer :: node

        base_bar = 0.0_dp
        enrichment_bar = 0.0_dp
        call validate_blending_inputs( &
            base_values, enriched_mask, enrichment_values, corrected_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(base_bar, 1) /= size(base_values, 1) .or. &
            size(base_bar, 2) /= size(base_values, 2) .or. &
            size(enrichment_bar) /= size(enrichment_values)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "corrected enrichment VJP has incompatible cotangents")
            return
        end if
        allocate(ramp(size(enrichment_values)))
        call blending_ramp(base_values, enriched_mask, ramp)
        do node = 1, size(base_values, 1)
            if (enriched_mask(node)) then
                base_bar(node, :) = corrected_bar*enrichment_values
            end if
        end do
        enrichment_bar = corrected_bar*ramp
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_blending_corrected_enrichment_vjp

    subroutine validate_blending_inputs( &
            base_values, enriched_mask, enrichment_values, corrected_values, status)
        real(dp), intent(in) :: base_values(:, :), enrichment_values(:)
        logical, intent(in) :: enriched_mask(:)
        real(dp), intent(in) :: corrected_values(:)
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "corrected enrichment has incompatible arrays")
        if (size(base_values, 1) < 1 .or. size(base_values, 2) < 1 .or. &
            size(enriched_mask) /= size(base_values, 1) .or. &
            size(enrichment_values) /= size(base_values, 2) .or. &
            size(corrected_values) /= size(base_values, 2)) return
        if (any(.not. ieee_is_finite(base_values)) .or. &
            any(.not. ieee_is_finite(enrichment_values)) .or. &
            any(.not. ieee_is_finite(corrected_values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "corrected enrichment received non-finite data")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_blending_inputs

    pure subroutine blending_ramp(base_values, enriched_mask, ramp)
        real(dp), intent(in) :: base_values(:, :)
        logical, intent(in) :: enriched_mask(:)
        real(dp), intent(out) :: ramp(:)
        integer :: node

        ramp = 0.0_dp
        do node = 1, size(base_values, 1)
            if (enriched_mask(node)) ramp = ramp + base_values(node, :)
        end do
    end subroutine blending_ramp

end module fortfem_xfem_blending_correction
