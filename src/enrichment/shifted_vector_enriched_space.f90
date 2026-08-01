module fortfem_shifted_vector_enriched_space
    !! Matrix-level shifted Heaviside enrichment for vector-valued bases.
    !!
    !! The first axis is the vector component, the second is the base-basis
    !! index, and the third is the quadrature or plotting point.  The scalar
    !! sign mask is shared by every component of one basis function:
    !!
    !!   N_enr(c,i,q) = N(c,i,q) [H(phi_q) - H(phi_i)].
    !!
    !! This is a composition layer before a caller applies covariant or
    !! contravariant Piola maps and an H(curl)/H(div) exact-sequence policy.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: evaluate_shifted_vector_enriched_space
    public :: evaluate_shifted_vector_enriched_space_jvp
    public :: evaluate_shifted_vector_enriched_space_vjp

contains

    subroutine evaluate_shifted_vector_enriched_space( &
            base_values, level_values, anchor_values, enriched_values, status)
        real(dp), intent(in) :: base_values(:, :, :), level_values(:), anchor_values(:)
        real(dp), intent(out) :: enriched_values(:, :, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: level_sign(:)
        integer :: component, basis

        enriched_values = 0.0_dp
        call validate_space_inputs( &
            base_values, level_values, anchor_values, enriched_values, status)
        if (status%code /= FORTSPARSE_OK) return
        allocate(level_sign(size(level_values)))
        level_sign = sign_indicator(level_values)
        do basis = 1, size(base_values, 2)
            do component = 1, size(base_values, 1)
                enriched_values(component, basis, :) = &
                    base_values(component, basis, :)* &
                    (level_sign - sign_value(anchor_values(basis)))
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_shifted_vector_enriched_space

    subroutine evaluate_shifted_vector_enriched_space_jvp( &
            base_values, level_values, anchor_values, base_dot, level_dot, &
            anchor_dot, enriched_dot, status)
        real(dp), intent(in) :: base_values(:, :, :), level_values(:), anchor_values(:)
        real(dp), intent(in) :: base_dot(:, :, :), level_dot(:), anchor_dot(:)
        real(dp), intent(out) :: enriched_dot(:, :, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: level_sign(:)
        integer :: component, basis

        enriched_dot = 0.0_dp
        call validate_space_inputs( &
            base_values, level_values, anchor_values, enriched_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (any(shape(base_dot) /= shape(base_values)) .or. &
            size(level_dot) /= size(level_values) .or. &
            size(anchor_dot) /= size(anchor_values) .or. &
            any(.not. ieee_is_finite(base_dot)) .or. &
            any(.not. ieee_is_finite(level_dot)) .or. &
            any(.not. ieee_is_finite(anchor_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "shifted vector enriched space JVP has incompatible increments")
            return
        end if
        allocate(level_sign(size(level_values)))
        level_sign = sign_indicator(level_values)
        do basis = 1, size(base_values, 2)
            do component = 1, size(base_values, 1)
                enriched_dot(component, basis, :) = &
                    base_dot(component, basis, :)* &
                    (level_sign - sign_value(anchor_values(basis)))
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_shifted_vector_enriched_space_jvp

    subroutine evaluate_shifted_vector_enriched_space_vjp( &
            base_values, level_values, anchor_values, enriched_bar, base_bar, &
            level_bar, anchor_bar, status)
        real(dp), intent(in) :: base_values(:, :, :), level_values(:), anchor_values(:)
        real(dp), intent(in) :: enriched_bar(:, :, :)
        real(dp), intent(out) :: base_bar(:, :, :), level_bar(:), anchor_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: level_sign(:)
        integer :: component, basis

        base_bar = 0.0_dp
        level_bar = 0.0_dp
        anchor_bar = 0.0_dp
        call validate_space_inputs( &
            base_values, level_values, anchor_values, enriched_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (any(shape(base_bar) /= shape(base_values)) .or. &
            size(level_bar) /= size(level_values) .or. &
            size(anchor_bar) /= size(anchor_values) .or. &
            any(.not. ieee_is_finite(enriched_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "shifted vector enriched space VJP has incompatible cotangents")
            return
        end if
        allocate(level_sign(size(level_values)))
        level_sign = sign_indicator(level_values)
        do basis = 1, size(base_values, 2)
            do component = 1, size(base_values, 1)
                base_bar(component, basis, :) = &
                    enriched_bar(component, basis, :)* &
                    (level_sign - sign_value(anchor_values(basis)))
            end do
        end do
        ! Fixed level-set and anchor signs carry no derivative.
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_shifted_vector_enriched_space_vjp

    subroutine validate_space_inputs( &
            base_values, level_values, anchor_values, enriched_values, status)
        real(dp), intent(in) :: base_values(:, :, :), level_values(:), anchor_values(:)
        real(dp), intent(in) :: enriched_values(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "shifted vector enriched space has incompatible arrays")
        if (size(base_values, 1) < 1 .or. size(base_values, 2) < 1 .or. &
            size(base_values, 3) < 1 .or. &
            size(level_values) /= size(base_values, 3) .or. &
            size(anchor_values) /= size(base_values, 2) .or. &
            any(shape(enriched_values) /= shape(base_values))) return
        if (any(.not. ieee_is_finite(base_values)) .or. &
            any(.not. ieee_is_finite(level_values)) .or. &
            any(.not. ieee_is_finite(anchor_values)) .or. &
            any(.not. ieee_is_finite(enriched_values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "shifted vector enriched space received non-finite data")
            return
        end if
        if (any(level_values == 0.0_dp) .or. &
            any(anchor_values == 0.0_dp)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "shifted vector enriched space encountered a topology event")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_space_inputs

    pure function sign_indicator(values) result(signs)
        real(dp), intent(in) :: values(:)
        real(dp) :: signs(size(values))
        integer :: item

        do item = 1, size(values)
            signs(item) = sign_value(values(item))
        end do
    end function sign_indicator

    pure real(dp) function sign_value(value)
        real(dp), intent(in) :: value

        if (value > 0.0_dp) then
            sign_value = 1.0_dp
        else
            sign_value = 0.0_dp
        end if
    end function sign_value

end module fortfem_shifted_vector_enriched_space
