module fortfem_heaviside_enrichment
    !! Piecewise-smooth shifted Heaviside enrichment for XFEM/GFEM bases.
    !!
    !! The enrichment is H(phi(x)) - H(phi(x_anchor)) with H=1 on the plus
    !! side and H=0 on the minus side.  A zero level value is a topology event
    !! and is rejected rather than assigned a misleading derivative.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_shifted_heaviside_enrichment
    public :: evaluate_shifted_heaviside_enrichment_jvp
    public :: evaluate_shifted_heaviside_enrichment_vjp

contains

    subroutine evaluate_shifted_heaviside_enrichment( &
            level_values, anchor_values, values, status)
        real(dp), intent(in) :: level_values(:), anchor_values(:)
        real(dp), intent(out) :: values(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: item

        values = 0.0_dp
        call validate_enrichment_inputs( &
            level_values, anchor_values, values, status)
        if (status%code /= FORTSPARSE_OK) return
        do item = 1, size(level_values)
            values(item) = heaviside(level_values(item)) - &
                heaviside(anchor_values(item))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_shifted_heaviside_enrichment

    subroutine evaluate_shifted_heaviside_enrichment_jvp( &
            level_values, anchor_values, level_dot, anchor_dot, values_dot, status)
        real(dp), intent(in) :: level_values(:), anchor_values(:)
        real(dp), intent(in) :: level_dot(:), anchor_dot(:)
        real(dp), intent(out) :: values_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        values_dot = 0.0_dp
        call validate_enrichment_inputs( &
            level_values, anchor_values, values_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(level_dot) /= size(level_values) .or. &
            size(anchor_dot) /= size(level_values) .or. &
            any(.not. ieee_is_finite(level_dot)) .or. &
            any(.not. ieee_is_finite(anchor_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "shifted Heaviside JVP has incompatible increments")
            return
        end if
        ! The activation signs are fixed; H is locally constant off phi=0.
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_shifted_heaviside_enrichment_jvp

    subroutine evaluate_shifted_heaviside_enrichment_vjp( &
            level_values, anchor_values, values_bar, level_bar, anchor_bar, status)
        real(dp), intent(in) :: level_values(:), anchor_values(:), values_bar(:)
        real(dp), intent(out) :: level_bar(:), anchor_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        level_bar = 0.0_dp
        anchor_bar = 0.0_dp
        call validate_enrichment_inputs( &
            level_values, anchor_values, values_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(level_bar) /= size(level_values) .or. &
            size(anchor_bar) /= size(level_values) .or. &
            any(.not. ieee_is_finite(values_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "shifted Heaviside VJP has incompatible cotangents")
            return
        end if
        ! The fixed-sign derivative is zero; topology events are rejected.
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_shifted_heaviside_enrichment_vjp

    subroutine validate_enrichment_inputs( &
            level_values, anchor_values, values, status)
        real(dp), intent(in) :: level_values(:), anchor_values(:), values(:)
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "shifted Heaviside enrichment has incompatible arrays")
        if (size(level_values) < 1 .or. &
            size(anchor_values) /= size(level_values) .or. &
            size(values) /= size(level_values)) return
        if (any(.not. ieee_is_finite(level_values)) .or. &
            any(.not. ieee_is_finite(anchor_values)) .or. &
            any(.not. ieee_is_finite(values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "shifted Heaviside enrichment received non-finite data")
            return
        end if
        if (any(level_values == 0.0_dp) .or. &
            any(anchor_values == 0.0_dp)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "shifted Heaviside enrichment encountered a topology event")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_enrichment_inputs

    pure real(dp) function heaviside(value)
        real(dp), intent(in) :: value

        if (value > 0.0_dp) then
            heaviside = 1.0_dp
        else
            heaviside = 0.0_dp
        end if
    end function heaviside

end module fortfem_heaviside_enrichment
