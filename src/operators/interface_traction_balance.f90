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
    public :: assemble_traction_jump
    public :: assemble_traction_jump_jvp
    public :: assemble_traction_jump_vjp

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

    subroutine assemble_traction_jump( &
            plus_traction, minus_traction, target, residual, status)
        !! Return the full vector traction jump minus a target traction.
        real(dp), intent(in) :: plus_traction(:), minus_traction(:), target(:)
        real(dp), intent(out) :: residual(:)
        type(fortsparse_status_t), intent(out) :: status

        residual = 0.0_dp
        if (.not. validate_vector_tractions( &
            plus_traction, minus_traction, target, residual, status)) return
        residual = plus_traction - minus_traction - target
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_traction_jump

    subroutine assemble_traction_jump_jvp( &
            plus_traction, minus_traction, target, plus_dot, minus_dot, &
            target_dot, residual_dot, status)
        !! Apply the product-rule JVP of the full traction jump.
        real(dp), intent(in) :: plus_traction(:), minus_traction(:), target(:)
        real(dp), intent(in) :: plus_dot(:), minus_dot(:), target_dot(:)
        real(dp), intent(out) :: residual_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        residual_dot = 0.0_dp
        if (.not. validate_vector_tractions( &
            plus_traction, minus_traction, target, residual_dot, status)) return
        if (size(plus_dot) /= 3 .or. size(minus_dot) /= 3 .or. &
            size(target_dot) /= 3 .or. any(.not. ieee_is_finite(plus_dot)) .or. &
            any(.not. ieee_is_finite(minus_dot)) .or. &
            any(.not. ieee_is_finite(target_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "vector traction jump JVP has incompatible increments")
            return
        end if
        residual_dot = plus_dot - minus_dot - target_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_traction_jump_jvp

    subroutine assemble_traction_jump_vjp( &
            plus_traction, minus_traction, target, residual_bar, plus_bar, &
            minus_bar, target_bar, status)
        !! Apply the real VJP of the full traction jump.
        real(dp), intent(in) :: plus_traction(:), minus_traction(:), target(:)
        real(dp), intent(in) :: residual_bar(:)
        real(dp), intent(out) :: plus_bar(:), minus_bar(:), target_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        plus_bar = 0.0_dp
        minus_bar = 0.0_dp
        target_bar = 0.0_dp
        if (.not. validate_vector_tractions( &
            plus_traction, minus_traction, target, residual_bar, status)) return
        if (size(plus_bar) /= 3 .or. size(minus_bar) /= 3 .or. &
            size(target_bar) /= 3 .or. any(.not. ieee_is_finite(residual_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "vector traction jump VJP has incompatible cotangents")
            return
        end if
        plus_bar = residual_bar
        minus_bar = -residual_bar
        target_bar = -residual_bar
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_traction_jump_vjp

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

    logical function validate_vector_tractions( &
            plus_traction, minus_traction, target, result, status) result(valid)
        real(dp), intent(in) :: plus_traction(:), minus_traction(:), target(:)
        real(dp), intent(in) :: result(:)
        type(fortsparse_status_t), intent(out) :: status

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "vector traction jump has incompatible arrays")
        if (size(plus_traction) /= 3 .or. size(minus_traction) /= 3 .or. &
            size(target) /= 3 .or. size(result) /= 3) return
        if (any(.not. ieee_is_finite(plus_traction)) .or. &
            any(.not. ieee_is_finite(minus_traction)) .or. &
            any(.not. ieee_is_finite(target)) .or. &
            any(.not. ieee_is_finite(result))) return
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_vector_tractions

end module fortfem_interface_traction_balance
