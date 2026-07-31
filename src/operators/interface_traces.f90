module fortfem_interface_traces
    !! Orientation-safe scalar and vector interface trace algebra.
    !!
    !! The supplied normal points from the minus side to the plus side.  All
    !! jumps therefore use plus-minus ordering.  Tangential traces are
    !! orthogonal projections, while the rotated jump n x (plus-minus) is the
    !! work-conjugate surface-current trace for Ampere interfaces.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    real(dp), parameter :: unit_tolerance = 1.0e-12_dp

    public :: compute_interface_scalar_jump_average
    public :: compute_interface_vector_traces

contains

    pure subroutine compute_interface_scalar_jump_average( &
            plus_value, minus_value, average, jump, status)
        real(dp), intent(in) :: plus_value, minus_value
        real(dp), intent(out) :: average, jump
        type(fortsparse_status_t), intent(out) :: status

        average = 0.0_dp
        jump = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "interface scalar trace received a non-finite value")
        if (.not. ieee_is_finite(plus_value) .or. &
            .not. ieee_is_finite(minus_value)) return
        average = 0.5_dp*(plus_value + minus_value)
        jump = plus_value - minus_value
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compute_interface_scalar_jump_average

    pure subroutine compute_interface_vector_traces( &
            plus_vector, minus_vector, normal, plus_normal, minus_normal, &
            normal_average, normal_jump, plus_tangent, minus_tangent, &
            tangent_average, tangent_jump, rotated_jump, status)
        real(dp), intent(in) :: plus_vector(:), minus_vector(:), normal(:)
        real(dp), intent(out) :: plus_normal, minus_normal
        real(dp), intent(out) :: normal_average, normal_jump
        real(dp), intent(out) :: plus_tangent(:), minus_tangent(:)
        real(dp), intent(out) :: tangent_average(:), tangent_jump(:)
        real(dp), intent(out) :: rotated_jump(:)
        type(fortsparse_status_t), intent(out) :: status

        plus_normal = 0.0_dp
        minus_normal = 0.0_dp
        normal_average = 0.0_dp
        normal_jump = 0.0_dp
        plus_tangent = 0.0_dp
        minus_tangent = 0.0_dp
        tangent_average = 0.0_dp
        tangent_jump = 0.0_dp
        rotated_jump = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "interface vector trace received incompatible arrays")
        if (size(plus_vector) /= 3 .or. size(minus_vector) /= 3 .or. &
            size(normal) /= 3 .or. size(plus_tangent) /= 3 .or. &
            size(minus_tangent) /= 3 .or. size(tangent_average) /= 3 .or. &
            size(tangent_jump) /= 3 .or. size(rotated_jump) /= 3) return
        if (any(.not. ieee_is_finite(plus_vector)) .or. &
            any(.not. ieee_is_finite(minus_vector)) .or. &
            any(.not. ieee_is_finite(normal))) return
        if (abs(dot_product(normal, normal) - 1.0_dp) > unit_tolerance) return

        plus_normal = dot_product(normal, plus_vector)
        minus_normal = dot_product(normal, minus_vector)
        normal_average = 0.5_dp*(plus_normal + minus_normal)
        normal_jump = plus_normal - minus_normal
        plus_tangent = plus_vector - plus_normal*normal
        minus_tangent = minus_vector - minus_normal*normal
        tangent_average = 0.5_dp*(plus_tangent + minus_tangent)
        tangent_jump = plus_tangent - minus_tangent
        call cross_product(normal, plus_vector - minus_vector, rotated_jump)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compute_interface_vector_traces

    pure subroutine cross_product(first, second, result)
        real(dp), intent(in) :: first(:), second(:)
        real(dp), intent(out) :: result(:)

        result(1) = first(2)*second(3) - first(3)*second(2)
        result(2) = first(3)*second(1) - first(1)*second(3)
        result(3) = first(1)*second(2) - first(2)*second(1)
    end subroutine cross_product

end module fortfem_interface_traces
