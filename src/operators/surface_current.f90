module fortfem_surface_current
    !! Orientation-safe surface current from an Ampere trace jump.
    !!
    !! The supplied normal points from the minus side to the plus side.  The
    !! surface current is K = n x (H_plus - H_minus), and the integrated
    !! ledger is the positive surface-quadrature sum of K.  This is a generic
    !! trace operation; constitutive and application boundary laws stay above
    !! FortFEM.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    real(dp), parameter :: unit_tolerance = 1.0e-12_dp

    public :: assemble_interface_surface_current
    public :: assemble_interface_surface_current_jvp
    public :: assemble_interface_surface_current_vjp

contains

    subroutine assemble_interface_surface_current( &
            plus_field, minus_field, normals, surface_weights, &
            surface_current, integrated_current, status)
        real(dp), intent(in) :: plus_field(:, :), minus_field(:, :)
        real(dp), intent(in) :: normals(:, :), surface_weights(:)
        real(dp), intent(out) :: surface_current(:, :), integrated_current(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: quadrature, quadrature_count
        real(dp) :: jump(3)

        surface_current = 0.0_dp
        integrated_current = 0.0_dp
        call validate_surface_current_inputs( &
            plus_field, minus_field, normals, surface_weights, &
            quadrature_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(surface_current, 1) /= quadrature_count .or. &
            size(surface_current, 2) /= 3 .or. &
            size(integrated_current) /= 3) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface current outputs have incompatible arrays")
            return
        end if

        do quadrature = 1, quadrature_count
            jump = plus_field(quadrature, :) - minus_field(quadrature, :)
            call cross_product(normals(quadrature, :), jump, &
                surface_current(quadrature, :))
            integrated_current = integrated_current + &
                surface_weights(quadrature)*surface_current(quadrature, :)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_interface_surface_current

    subroutine assemble_interface_surface_current_jvp( &
            plus_field, minus_field, normals, surface_weights, plus_dot, &
            minus_dot, normals_dot, surface_weights_dot, surface_current_dot, &
            integrated_current_dot, status)
        real(dp), intent(in) :: plus_field(:, :), minus_field(:, :)
        real(dp), intent(in) :: normals(:, :), surface_weights(:)
        real(dp), intent(in) :: plus_dot(:, :), minus_dot(:, :), normals_dot(:, :)
        real(dp), intent(in) :: surface_weights_dot(:)
        real(dp), intent(out) :: surface_current_dot(:, :)
        real(dp), intent(out) :: integrated_current_dot(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: quadrature, quadrature_count
        real(dp) :: jump(3), jump_dot(3), current(3), current_dot(3)
        real(dp) :: normal_dot_term(3), jump_dot_term(3)

        surface_current_dot = 0.0_dp
        integrated_current_dot = 0.0_dp
        call validate_surface_current_inputs( &
            plus_field, minus_field, normals, surface_weights, &
            quadrature_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(plus_dot, 1) /= quadrature_count .or. &
            size(plus_dot, 2) /= 3 .or. &
            size(minus_dot, 1) /= quadrature_count .or. &
            size(minus_dot, 2) /= 3 .or. &
            size(normals_dot, 1) /= quadrature_count .or. &
            size(normals_dot, 2) /= 3 .or. &
            size(surface_weights_dot) /= quadrature_count .or. &
            size(surface_current_dot, 1) /= quadrature_count .or. &
            size(surface_current_dot, 2) /= 3 .or. &
            size(integrated_current_dot) /= 3) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface current JVP received incompatible arrays")
            return
        end if
        if (any(.not. ieee_is_finite(plus_dot)) .or. &
            any(.not. ieee_is_finite(minus_dot)) .or. &
            any(.not. ieee_is_finite(normals_dot)) .or. &
            any(.not. ieee_is_finite(surface_weights_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface current JVP received non-finite increments")
            return
        end if

        do quadrature = 1, quadrature_count
            jump = plus_field(quadrature, :) - minus_field(quadrature, :)
            jump_dot = plus_dot(quadrature, :) - minus_dot(quadrature, :)
            call cross_product(normals_dot(quadrature, :), jump, &
                normal_dot_term)
            call cross_product(normals(quadrature, :), jump_dot, &
                jump_dot_term)
            current_dot = normal_dot_term + jump_dot_term
            surface_current_dot(quadrature, :) = current_dot
            call cross_product(normals(quadrature, :), jump, current)
            integrated_current_dot = integrated_current_dot + &
                surface_weights_dot(quadrature)*current + &
                surface_weights(quadrature)*current_dot
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_interface_surface_current_jvp

    subroutine assemble_interface_surface_current_vjp( &
            plus_field, minus_field, normals, surface_weights, surface_current_bar, &
            integrated_current_bar, plus_bar, minus_bar, normals_bar, &
            surface_weights_bar, status)
        real(dp), intent(in) :: plus_field(:, :), minus_field(:, :)
        real(dp), intent(in) :: normals(:, :), surface_weights(:)
        real(dp), intent(in) :: surface_current_bar(:, :)
        real(dp), intent(in) :: integrated_current_bar(:)
        real(dp), intent(out) :: plus_bar(:, :), minus_bar(:, :), normals_bar(:, :)
        real(dp), intent(out) :: surface_weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: quadrature, quadrature_count
        real(dp) :: jump(3), current(3), current_bar_total(3)

        plus_bar = 0.0_dp
        minus_bar = 0.0_dp
        normals_bar = 0.0_dp
        surface_weights_bar = 0.0_dp
        call validate_surface_current_inputs( &
            plus_field, minus_field, normals, surface_weights, &
            quadrature_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(surface_current_bar, 1) /= quadrature_count .or. &
            size(surface_current_bar, 2) /= 3 .or. &
            size(integrated_current_bar) /= 3 .or. &
            size(plus_bar, 1) /= quadrature_count .or. size(plus_bar, 2) /= 3 .or. &
            size(minus_bar, 1) /= quadrature_count .or. &
            size(minus_bar, 2) /= 3 .or. &
            size(normals_bar, 1) /= quadrature_count .or. &
            size(normals_bar, 2) /= 3 .or. &
            size(surface_weights_bar) /= quadrature_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface current VJP received incompatible arrays")
            return
        end if
        if (any(.not. ieee_is_finite(surface_current_bar)) .or. &
            any(.not. ieee_is_finite(integrated_current_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface current VJP received non-finite cotangents")
            return
        end if

        do quadrature = 1, quadrature_count
            jump = plus_field(quadrature, :) - minus_field(quadrature, :)
            call cross_product(normals(quadrature, :), jump, current)
            current_bar_total = surface_current_bar(quadrature, :) + &
                surface_weights(quadrature)*integrated_current_bar
            call cross_product(current_bar_total, normals(quadrature, :), &
                plus_bar(quadrature, :))
            minus_bar(quadrature, :) = -plus_bar(quadrature, :)
            call cross_product(jump, current_bar_total, &
                normals_bar(quadrature, :))
            surface_weights_bar(quadrature) = dot_product( &
                integrated_current_bar, current)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_interface_surface_current_vjp

    subroutine validate_surface_current_inputs( &
            plus_field, minus_field, normals, surface_weights, &
            quadrature_count, status)
        real(dp), intent(in) :: plus_field(:, :), minus_field(:, :)
        real(dp), intent(in) :: normals(:, :), surface_weights(:)
        integer, intent(out) :: quadrature_count
        type(fortsparse_status_t), intent(out) :: status
        integer :: quadrature

        quadrature_count = size(plus_field, 1)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "surface current received incompatible arrays")
        if (quadrature_count < 1) return
        if (size(plus_field, 2) /= 3 .or. &
            size(minus_field, 1) /= quadrature_count .or. &
            size(minus_field, 2) /= 3 .or. &
            size(normals, 1) /= quadrature_count .or. &
            size(normals, 2) /= 3 .or. &
            size(surface_weights) /= quadrature_count) return
        if (any(.not. ieee_is_finite(plus_field)) .or. &
            any(.not. ieee_is_finite(minus_field)) .or. &
            any(.not. ieee_is_finite(normals)) .or. &
            any(.not. ieee_is_finite(surface_weights))) return
        if (any(surface_weights <= 0.0_dp)) return
        do quadrature = 1, quadrature_count
            if (abs(dot_product(normals(quadrature, :), &
                normals(quadrature, :)) - 1.0_dp) > unit_tolerance) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_surface_current_inputs

    pure subroutine cross_product(first, second, result)
        real(dp), intent(in) :: first(:), second(:)
        real(dp), intent(out) :: result(:)

        result(1) = first(2)*second(3) - first(3)*second(2)
        result(2) = first(3)*second(1) - first(1)*second(3)
        result(3) = first(1)*second(2) - first(2)*second(1)
    end subroutine cross_product

end module fortfem_surface_current
