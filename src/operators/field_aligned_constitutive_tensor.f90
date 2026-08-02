module fortfem_field_aligned_constitutive_tensor
    !! Pointwise field-aligned symmetric and gyrotropic constitutive tensor.
    !!
    !! For a unit direction b, the tensor is
    !!
    !!   K = k_perpendicular I + (k_parallel-k_perpendicular) b b^T
    !!       + k_hall [b]_x,
    !!
    !! where [b]_x v = b x v.  The Hall coefficient is an optional trailing
    !! argument and defaults to zero.  This is a caller-owned constitutive
    !! ingredient; it does not select a plasma closure or a material model.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_cgl_pressure_tensor, only: &
        evaluate_cgl_pressure_tensor, evaluate_cgl_pressure_tensor_jvp, &
        evaluate_cgl_pressure_tensor_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    real(dp), parameter :: unit_tolerance = 1.0e-12_dp

    public :: evaluate_field_aligned_constitutive_tensor
    public :: evaluate_field_aligned_constitutive_tensor_jvp
    public :: evaluate_field_aligned_constitutive_tensor_vjp

contains

    pure subroutine evaluate_field_aligned_constitutive_tensor( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            tensor, status, hall_coefficient)
        real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient
        real(dp), intent(in) :: unit_direction(:)
        real(dp), intent(out) :: tensor(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: hall_coefficient

        real(dp) :: hall, symmetric(3, 3), cross(3, 3)

        tensor = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "field-aligned constitutive tensor received invalid arrays")
        hall = 0.0_dp
        if (present(hall_coefficient)) hall = hall_coefficient
        if (.not. valid_inputs(parallel_coefficient, perpendicular_coefficient, &
            hall, unit_direction, tensor)) return
        call evaluate_cgl_pressure_tensor( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            symmetric, status)
        if (status%code /= FORTSPARSE_OK) return
        call cross_product_matrix(unit_direction, cross)
        tensor = symmetric + hall*cross
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_field_aligned_constitutive_tensor

    pure subroutine evaluate_field_aligned_constitutive_tensor_jvp( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            parallel_coefficient_dot, perpendicular_coefficient_dot, &
            direction_dot, tensor_dot, status, hall_coefficient, &
            hall_coefficient_dot)
        real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient
        real(dp), intent(in) :: unit_direction(:)
        real(dp), intent(in) :: parallel_coefficient_dot
        real(dp), intent(in) :: perpendicular_coefficient_dot
        real(dp), intent(in) :: direction_dot(:)
        real(dp), intent(out) :: tensor_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: hall_coefficient, hall_coefficient_dot

        real(dp) :: hall, hall_dot, symmetric_dot(3, 3), cross(3, 3)
        real(dp) :: cross_dot(3, 3)

        tensor_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "field-aligned constitutive tensor JVP received invalid arrays")
        hall = 0.0_dp
        hall_dot = 0.0_dp
        if (present(hall_coefficient)) hall = hall_coefficient
        if (present(hall_coefficient_dot)) hall_dot = hall_coefficient_dot
        if (.not. valid_inputs(parallel_coefficient, perpendicular_coefficient, &
            hall, unit_direction, tensor_dot)) return
        if (size(direction_dot) /= 3 .or. &
            any(.not. ieee_is_finite(direction_dot)) .or. &
            .not. ieee_is_finite(parallel_coefficient_dot) .or. &
            .not. ieee_is_finite(perpendicular_coefficient_dot) .or. &
            .not. ieee_is_finite(hall_dot)) return
        call evaluate_cgl_pressure_tensor_jvp( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            parallel_coefficient_dot, perpendicular_coefficient_dot, &
            direction_dot, symmetric_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        call cross_product_matrix(unit_direction, cross)
        call cross_product_matrix(direction_dot, cross_dot)
        tensor_dot = symmetric_dot + hall_dot*cross + hall*cross_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_field_aligned_constitutive_tensor_jvp

    pure subroutine evaluate_field_aligned_constitutive_tensor_vjp( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            tensor_bar, parallel_coefficient_bar, perpendicular_coefficient_bar, &
            direction_bar, status, hall_coefficient, hall_coefficient_bar)
        real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient
        real(dp), intent(in) :: unit_direction(:), tensor_bar(:, :)
        real(dp), intent(out) :: parallel_coefficient_bar
        real(dp), intent(out) :: perpendicular_coefficient_bar
        real(dp), intent(out) :: direction_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: hall_coefficient
        real(dp), intent(out), optional :: hall_coefficient_bar

        real(dp) :: hall, symmetric_direction_bar(3), cross(3, 3)

        parallel_coefficient_bar = 0.0_dp
        perpendicular_coefficient_bar = 0.0_dp
        direction_bar = 0.0_dp
        if (present(hall_coefficient_bar)) hall_coefficient_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "field-aligned constitutive tensor VJP received invalid arrays")
        hall = 0.0_dp
        if (present(hall_coefficient)) hall = hall_coefficient
        if (.not. valid_inputs(parallel_coefficient, perpendicular_coefficient, &
            hall, unit_direction, tensor_bar)) return
        if (size(direction_bar) /= 3) return
        if (any(.not. ieee_is_finite(tensor_bar))) return

        call evaluate_cgl_pressure_tensor_vjp( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            tensor_bar, parallel_coefficient_bar, perpendicular_coefficient_bar, &
            symmetric_direction_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        direction_bar = symmetric_direction_bar
        call cross_product_matrix(unit_direction, cross)
        if (present(hall_coefficient_bar)) then
            hall_coefficient_bar = sum(tensor_bar*cross)
        end if
        direction_bar(1) = direction_bar(1) + hall*(tensor_bar(3, 2) - &
            tensor_bar(2, 3))
        direction_bar(2) = direction_bar(2) + hall*(tensor_bar(1, 3) - &
            tensor_bar(3, 1))
        direction_bar(3) = direction_bar(3) + hall*(tensor_bar(2, 1) - &
            tensor_bar(1, 2))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_field_aligned_constitutive_tensor_vjp

    pure logical function valid_inputs( &
            parallel_coefficient, perpendicular_coefficient, hall, &
            unit_direction, tensor) result(valid)
        real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient
        real(dp), intent(in) :: hall, unit_direction(:), tensor(:, :)

        valid = ieee_is_finite(parallel_coefficient) .and. &
            ieee_is_finite(perpendicular_coefficient) .and. &
            ieee_is_finite(hall) .and. size(unit_direction) == 3 .and. &
            size(tensor, 1) == 3 .and. size(tensor, 2) == 3 .and. &
            all(ieee_is_finite(unit_direction)) .and. &
            abs(dot_product(unit_direction, unit_direction) - 1.0_dp) <= &
            unit_tolerance
    end function valid_inputs

    pure subroutine cross_product_matrix(direction, matrix)
        real(dp), intent(in) :: direction(3)
        real(dp), intent(out) :: matrix(3, 3)

        matrix = 0.0_dp
        matrix(1, 2) = -direction(3)
        matrix(1, 3) = direction(2)
        matrix(2, 1) = direction(3)
        matrix(2, 3) = -direction(1)
        matrix(3, 1) = -direction(2)
        matrix(3, 2) = direction(1)
    end subroutine cross_product_matrix

end module fortfem_field_aligned_constitutive_tensor
