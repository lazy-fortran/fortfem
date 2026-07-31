module fortfem_shifted_enriched_basis
    !! Partition-of-unity composition of a base value and shifted activation.
    !!
    !! The basis value is N_i*(H(phi)-H(phi_i)).  This layer keeps the base
    !! approximation separate from cut-cell geometry and exposes the product
    !! rule needed by scalar XFEM/GFEM assembly.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_heaviside_enrichment, only: &
        evaluate_shifted_heaviside_enrichment, &
        evaluate_shifted_heaviside_enrichment_jvp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_shifted_enriched_basis
    public :: evaluate_shifted_enriched_basis_jvp
    public :: evaluate_shifted_enriched_basis_vjp

contains

    subroutine evaluate_shifted_enriched_basis( &
            base_values, level_values, anchor_values, enriched_values, status)
        real(dp), intent(in) :: base_values(:), level_values(:), anchor_values(:)
        real(dp), intent(out) :: enriched_values(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: activation(:)

        enriched_values = 0.0_dp
        call validate_basis_inputs( &
            base_values, level_values, anchor_values, enriched_values, status)
        if (status%code /= FORTSPARSE_OK) return
        allocate(activation(size(base_values)))
        call evaluate_shifted_heaviside_enrichment( &
            level_values, anchor_values, activation, status)
        if (status%code /= FORTSPARSE_OK) return
        enriched_values = base_values*activation
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_shifted_enriched_basis

    subroutine evaluate_shifted_enriched_basis_jvp( &
            base_values, level_values, anchor_values, base_dot, level_dot, &
            anchor_dot, enriched_dot, status)
        real(dp), intent(in) :: base_values(:), level_values(:), anchor_values(:)
        real(dp), intent(in) :: base_dot(:), level_dot(:), anchor_dot(:)
        real(dp), intent(out) :: enriched_dot(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: activation(:), activation_dot(:)

        enriched_dot = 0.0_dp
        call validate_basis_inputs( &
            base_values, level_values, anchor_values, enriched_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(base_dot) /= size(base_values) .or. &
            size(level_dot) /= size(base_values) .or. &
            size(anchor_dot) /= size(base_values) .or. &
            any(.not. ieee_is_finite(base_dot)) .or. &
            any(.not. ieee_is_finite(level_dot)) .or. &
            any(.not. ieee_is_finite(anchor_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "shifted enriched basis JVP has incompatible increments")
            return
        end if
        allocate(activation(size(base_values)), activation_dot(size(base_values)))
        call evaluate_shifted_heaviside_enrichment( &
            level_values, anchor_values, activation, status)
        if (status%code /= FORTSPARSE_OK) return
        call evaluate_shifted_heaviside_enrichment_jvp( &
            level_values, anchor_values, level_dot, anchor_dot, activation_dot, &
            status)
        if (status%code /= FORTSPARSE_OK) return
        enriched_dot = base_dot*activation + base_values*activation_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_shifted_enriched_basis_jvp

    subroutine evaluate_shifted_enriched_basis_vjp( &
            base_values, level_values, anchor_values, enriched_bar, base_bar, &
            level_bar, anchor_bar, status)
        real(dp), intent(in) :: base_values(:), level_values(:), anchor_values(:)
        real(dp), intent(in) :: enriched_bar(:)
        real(dp), intent(out) :: base_bar(:), level_bar(:), anchor_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: activation(:)

        base_bar = 0.0_dp
        level_bar = 0.0_dp
        anchor_bar = 0.0_dp
        call validate_basis_inputs( &
            base_values, level_values, anchor_values, enriched_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(base_bar) /= size(base_values) .or. &
            size(level_bar) /= size(base_values) .or. &
            size(anchor_bar) /= size(base_values) .or. &
            any(.not. ieee_is_finite(enriched_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "shifted enriched basis VJP has incompatible cotangents")
            return
        end if
        allocate(activation(size(base_values)))
        call evaluate_shifted_heaviside_enrichment( &
            level_values, anchor_values, activation, status)
        if (status%code /= FORTSPARSE_OK) return
        base_bar = enriched_bar*activation
        ! Fixed activation signs carry no derivative; topology events reject.
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_shifted_enriched_basis_vjp

    subroutine validate_basis_inputs( &
            base_values, level_values, anchor_values, enriched_values, status)
        real(dp), intent(in) :: base_values(:), level_values(:), anchor_values(:)
        real(dp), intent(in) :: enriched_values(:)
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "shifted enriched basis has incompatible arrays")
        if (size(base_values) < 1 .or. &
            size(level_values) /= size(base_values) .or. &
            size(anchor_values) /= size(base_values) .or. &
            size(enriched_values) /= size(base_values)) return
        if (any(.not. ieee_is_finite(base_values)) .or. &
            any(.not. ieee_is_finite(level_values)) .or. &
            any(.not. ieee_is_finite(anchor_values)) .or. &
            any(.not. ieee_is_finite(enriched_values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "shifted enriched basis received non-finite data")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_basis_inputs

end module fortfem_shifted_enriched_basis
