module fortfem_interface_traction_balance
    !! Neutral normal-traction jump residual for fitted or cut interfaces.
    !!
    !! The caller supplies the two physical traction vectors.  This block
    !! projects their jump onto an oriented unit normal and subtracts a
    !! caller-owned scalar target, so tensor-valued pressure, elasticity, or
    !! Maxwell-stress blocks can share one differentiable interface contract.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    real(dp), parameter :: unit_tolerance = 1.0e-12_dp

    public :: assemble_normal_traction_jump
    public :: assemble_normal_traction_jump_jvp
    public :: assemble_normal_traction_jump_vjp

contains

    subroutine assemble_normal_traction_jump( &
            plus_traction, minus_traction, normal, target, residual, status)
        real(dp), intent(in) :: plus_traction(:), minus_traction(:), normal(:)
        real(dp), intent(in) :: target
        real(dp), intent(out) :: residual
        type(fortsparse_status_t), intent(out) :: status

        residual = 0.0_dp
        call validate_traction_inputs( &
            plus_traction, minus_traction, normal, target, status)
        if (status%code /= FORTSPARSE_OK) return
        residual = dot_product(normal, plus_traction - minus_traction) - target
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_normal_traction_jump

    subroutine assemble_normal_traction_jump_jvp( &
            plus_traction, minus_traction, normal, target, plus_dot, minus_dot, &
            normal_dot, target_dot, residual_dot, status)
        real(dp), intent(in) :: plus_traction(:), minus_traction(:), normal(:)
        real(dp), intent(in) :: target, plus_dot(:), minus_dot(:), normal_dot(:)
        real(dp), intent(in) :: target_dot
        real(dp), intent(out) :: residual_dot
        type(fortsparse_status_t), intent(out) :: status

        residual_dot = 0.0_dp
        call validate_traction_inputs( &
            plus_traction, minus_traction, normal, target, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(plus_dot) /= 3 .or. size(minus_dot) /= 3 .or. &
            size(normal_dot) /= 3 .or. any(.not. ieee_is_finite(plus_dot)) .or. &
            any(.not. ieee_is_finite(minus_dot)) .or. &
            any(.not. ieee_is_finite(normal_dot)) .or. &
            .not. ieee_is_finite(target_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "normal traction jump JVP has incompatible increments")
            return
        end if
        residual_dot = dot_product(normal_dot, plus_traction - minus_traction) + &
            dot_product(normal, plus_dot - minus_dot) - target_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_normal_traction_jump_jvp

    subroutine assemble_normal_traction_jump_vjp( &
            plus_traction, minus_traction, normal, target, residual_bar, &
            plus_bar, minus_bar, normal_bar, target_bar, status)
        real(dp), intent(in) :: plus_traction(:), minus_traction(:), normal(:)
        real(dp), intent(in) :: target, residual_bar
        real(dp), intent(out) :: plus_bar(:), minus_bar(:), normal_bar(:)
        real(dp), intent(out) :: target_bar
        type(fortsparse_status_t), intent(out) :: status

        plus_bar = 0.0_dp
        minus_bar = 0.0_dp
        normal_bar = 0.0_dp
        target_bar = 0.0_dp
        call validate_traction_inputs( &
            plus_traction, minus_traction, normal, target, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(plus_bar) /= 3 .or. size(minus_bar) /= 3 .or. &
            size(normal_bar) /= 3 .or. .not. ieee_is_finite(residual_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "normal traction jump VJP has incompatible cotangents")
            return
        end if
        plus_bar = residual_bar*normal
        minus_bar = -plus_bar
        normal_bar = residual_bar*(plus_traction - minus_traction)
        target_bar = -residual_bar
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_normal_traction_jump_vjp

    subroutine validate_traction_inputs( &
            plus_traction, minus_traction, normal, target, status)
        real(dp), intent(in) :: plus_traction(:), minus_traction(:), normal(:)
        real(dp), intent(in) :: target
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "normal traction jump has incompatible arrays")
        if (size(plus_traction) /= 3 .or. size(minus_traction) /= 3 .or. &
            size(normal) /= 3) return
        if (any(.not. ieee_is_finite(plus_traction)) .or. &
            any(.not. ieee_is_finite(minus_traction)) .or. &
            any(.not. ieee_is_finite(normal)) .or. &
            .not. ieee_is_finite(target)) return
        if (abs(dot_product(normal, normal) - 1.0_dp) > unit_tolerance) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_traction_inputs

end module fortfem_interface_traction_balance
