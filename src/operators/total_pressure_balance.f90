module fortfem_total_pressure_balance
    !! Neutral total-pressure jump block for interface compositions.
    !!
    !! The caller owns the scalar pressure and magnetic-field values on both
    !! sides.  This module only evaluates the algebraic jump
    !! p_plus + |B_plus|^2/(2 mu) - p_minus - |B_minus|^2/(2 mu) - target.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_total_pressure_jump
    public :: assemble_total_pressure_jump_jvp
    public :: assemble_total_pressure_jump_vjp

contains

    subroutine assemble_total_pressure_jump( &
            plus_pressure, minus_pressure, plus_field, minus_field, permeability, &
            target, residual, status)
        real(dp), intent(in) :: plus_pressure, minus_pressure
        real(dp), intent(in) :: plus_field(:), minus_field(:), permeability, target
        real(dp), intent(out) :: residual
        type(fortsparse_status_t), intent(out) :: status

        residual = 0.0_dp
        if (.not. valid_inputs(plus_pressure, minus_pressure, plus_field, &
            minus_field, permeability, target)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "total-pressure jump received invalid inputs")
            return
        end if
        residual = plus_pressure + 0.5_dp*dot_product(plus_field, plus_field)/ &
            permeability - minus_pressure - &
            0.5_dp*dot_product(minus_field, minus_field)/permeability - target
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_total_pressure_jump

    subroutine assemble_total_pressure_jump_jvp( &
            plus_pressure, minus_pressure, plus_field, minus_field, permeability, &
            target, plus_pressure_dot, minus_pressure_dot, plus_field_dot, &
            minus_field_dot, permeability_dot, target_dot, residual_dot, status)
        real(dp), intent(in) :: plus_pressure, minus_pressure
        real(dp), intent(in) :: plus_field(:), minus_field(:), permeability, target
        real(dp), intent(in) :: plus_pressure_dot, minus_pressure_dot
        real(dp), intent(in) :: plus_field_dot(:), minus_field_dot(:)
        real(dp), intent(in) :: permeability_dot, target_dot
        real(dp), intent(out) :: residual_dot
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: residual

        residual_dot = 0.0_dp
        call assemble_total_pressure_jump(plus_pressure, minus_pressure, plus_field, &
            minus_field, permeability, target, residual, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(plus_field_dot) /= 3 .or. size(minus_field_dot) /= 3 .or. &
            .not. all(ieee_is_finite([plus_pressure_dot, minus_pressure_dot, &
            permeability_dot, target_dot])) .or. &
            .not. all(ieee_is_finite(plus_field_dot)) .or. &
            .not. all(ieee_is_finite(minus_field_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "total-pressure jump JVP has invalid increments")
            return
        end if
        residual_dot = plus_pressure_dot - minus_pressure_dot + &
            dot_product(plus_field, plus_field_dot)/permeability - &
            dot_product(minus_field, minus_field_dot)/permeability - &
            0.5_dp*(dot_product(plus_field, plus_field) - &
            dot_product(minus_field, minus_field))*permeability_dot/ &
            permeability**2 - target_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_total_pressure_jump_jvp

    subroutine assemble_total_pressure_jump_vjp( &
            plus_pressure, minus_pressure, plus_field, minus_field, permeability, &
            target, residual_bar, plus_pressure_bar, minus_pressure_bar, &
            plus_field_bar, minus_field_bar, permeability_bar, target_bar, status)
        real(dp), intent(in) :: plus_pressure, minus_pressure
        real(dp), intent(in) :: plus_field(:), minus_field(:), permeability, target
        real(dp), intent(in) :: residual_bar
        real(dp), intent(out) :: plus_pressure_bar, minus_pressure_bar
        real(dp), intent(out) :: plus_field_bar(:), minus_field_bar(:)
        real(dp), intent(out) :: permeability_bar, target_bar
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: residual

        plus_pressure_bar = 0.0_dp
        minus_pressure_bar = 0.0_dp
        plus_field_bar = 0.0_dp
        minus_field_bar = 0.0_dp
        permeability_bar = 0.0_dp
        target_bar = 0.0_dp
        call assemble_total_pressure_jump(plus_pressure, minus_pressure, plus_field, &
            minus_field, permeability, target, residual, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. ieee_is_finite(residual_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "total-pressure jump VJP has invalid cotangent")
            return
        end if
        plus_pressure_bar = residual_bar
        minus_pressure_bar = -residual_bar
        plus_field_bar = residual_bar*plus_field/permeability
        minus_field_bar = -residual_bar*minus_field/permeability
        permeability_bar = -0.5_dp*residual_bar*(dot_product(plus_field, plus_field) - &
            dot_product(minus_field, minus_field))/permeability**2
        target_bar = -residual_bar
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_total_pressure_jump_vjp

    logical function valid_inputs(plus_pressure, minus_pressure, plus_field, &
            minus_field, permeability, target) result(valid)
        real(dp), intent(in) :: plus_pressure, minus_pressure
        real(dp), intent(in) :: plus_field(:), minus_field(:), permeability, target

        valid = size(plus_field) == 3 .and. size(minus_field) == 3 .and. &
            ieee_is_finite(plus_pressure) .and. ieee_is_finite(minus_pressure) .and. &
            ieee_is_finite(permeability) .and. permeability > 0.0_dp .and. &
            ieee_is_finite(target) .and. all(ieee_is_finite(plus_field)) .and. &
            all(ieee_is_finite(minus_field))
    end function valid_inputs

end module fortfem_total_pressure_balance
