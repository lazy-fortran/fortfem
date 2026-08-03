module fortfem_field_aligned_tensor_pullback
    !! Metric pullback for a physical field-aligned constitutive tensor.
    !!
    !! For an oriented map x = F(xi), J = dF/dxi, and a physical tensor K,
    !! the scalar H1 diffusion contraction is represented in reference
    !! coordinates by
    !!
    !!   K_ref = det(J) J^{-1} K J^{-T}.
    !!
    !! The constitutive tensor is supplied by the neutral field-aligned
    !! operator.  This module owns only metric transformation and its exact
    !! tangent/transpose products; callers own quadrature and field physics.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_field_aligned_constitutive_tensor, only: &
        evaluate_field_aligned_constitutive_tensor, &
        evaluate_field_aligned_constitutive_tensor_jvp, &
        evaluate_field_aligned_constitutive_tensor_vjp
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3, det3_jvp, det3_vjp, inv3, inv3_jvp, inv3_vjp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: pullback_field_aligned_tensor
    public :: pullback_field_aligned_tensor_jvp
    public :: pullback_field_aligned_tensor_vjp

contains

    pure subroutine pullback_field_aligned_tensor( &
            jacobian, parallel_coefficient, perpendicular_coefficient, &
            unit_direction, reference_tensor, status, hall_coefficient)
        real(dp), intent(in) :: jacobian(3, 3)
        real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient
        real(dp), intent(in) :: unit_direction(:)
        real(dp), intent(out) :: reference_tensor(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: hall_coefficient

        real(dp) :: physical_tensor(3, 3), inverse_jacobian(3, 3)
        real(dp) :: determinant
        integer :: inverse_status

        reference_tensor = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "field-aligned tensor pullback received invalid arrays")
        if (size(reference_tensor, 1) /= 3 .or. &
            size(reference_tensor, 2) /= 3 .or. .not. valid_jacobian(jacobian)) return
        call evaluate_constitutive( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            physical_tensor, status, hall_coefficient)
        if (status%code /= FORTSPARSE_OK) return
        determinant = det3(jacobian)
        call inv3(jacobian, inverse_jacobian, inverse_status)
        if (inverse_status /= 0) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "field-aligned tensor pullback map is singular")
            return
        end if
        reference_tensor = determinant*matmul( &
            inverse_jacobian, matmul(physical_tensor, transpose(inverse_jacobian)))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine pullback_field_aligned_tensor

    pure subroutine pullback_field_aligned_tensor_jvp( &
            jacobian, parallel_coefficient, perpendicular_coefficient, &
            unit_direction, jacobian_dot, parallel_coefficient_dot, &
            perpendicular_coefficient_dot, direction_dot, reference_tensor_dot, &
            status, hall_coefficient, hall_coefficient_dot)
        real(dp), intent(in) :: jacobian(3, 3)
        real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient
        real(dp), intent(in) :: unit_direction(:), jacobian_dot(3, 3)
        real(dp), intent(in) :: parallel_coefficient_dot
        real(dp), intent(in) :: perpendicular_coefficient_dot
        real(dp), intent(in) :: direction_dot(:)
        real(dp), intent(out) :: reference_tensor_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: hall_coefficient, hall_coefficient_dot

        real(dp) :: physical_tensor(3, 3), physical_tensor_dot(3, 3)
        real(dp) :: inverse_jacobian(3, 3), inverse_jacobian_dot(3, 3)
        real(dp) :: determinant, determinant_dot
        integer :: inverse_status

        reference_tensor_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "field-aligned tensor pullback JVP received invalid arrays")
        if (size(reference_tensor_dot, 1) /= 3 .or. &
            size(reference_tensor_dot, 2) /= 3 .or. .not. valid_jacobian(jacobian)) return
        if (.not. all(ieee_is_finite(jacobian_dot)) .or. size(direction_dot) /= 3) return
        if (.not. ieee_is_finite(parallel_coefficient_dot) .or. &
            .not. ieee_is_finite(perpendicular_coefficient_dot)) return
        call evaluate_constitutive( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            physical_tensor, status, hall_coefficient)
        if (status%code /= FORTSPARSE_OK) return
        call evaluate_constitutive_jvp( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            parallel_coefficient_dot, perpendicular_coefficient_dot, direction_dot, &
            physical_tensor_dot, status, hall_coefficient, hall_coefficient_dot)
        if (status%code /= FORTSPARSE_OK) return
        call inv3_jvp( &
            jacobian, jacobian_dot, inverse_jacobian, inverse_jacobian_dot, &
            inverse_status)
        if (inverse_status /= 0) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "field-aligned tensor pullback JVP map is singular")
            return
        end if
        determinant = det3(jacobian)
        call det3_jvp(jacobian, jacobian_dot, determinant_dot)
        reference_tensor_dot = determinant_dot*matmul( &
            inverse_jacobian, matmul(physical_tensor, transpose(inverse_jacobian))) + &
            determinant*( &
            matmul(inverse_jacobian_dot, matmul(physical_tensor, transpose(inverse_jacobian))) + &
            matmul(inverse_jacobian, matmul(physical_tensor_dot, transpose(inverse_jacobian))) + &
            matmul(inverse_jacobian, matmul(physical_tensor, transpose(inverse_jacobian_dot))))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine pullback_field_aligned_tensor_jvp

    pure subroutine pullback_field_aligned_tensor_vjp( &
            jacobian, parallel_coefficient, perpendicular_coefficient, &
            unit_direction, reference_tensor_bar, jacobian_bar, &
            parallel_coefficient_bar, perpendicular_coefficient_bar, direction_bar, &
            status, hall_coefficient, hall_coefficient_bar)
        real(dp), intent(in) :: jacobian(3, 3)
        real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient
        real(dp), intent(in) :: unit_direction(:), reference_tensor_bar(:, :)
        real(dp), intent(out) :: jacobian_bar(3, 3)
        real(dp), intent(out) :: parallel_coefficient_bar
        real(dp), intent(out) :: perpendicular_coefficient_bar
        real(dp), intent(out) :: direction_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: hall_coefficient
        real(dp), intent(out), optional :: hall_coefficient_bar

        real(dp) :: physical_tensor(3, 3), physical_tensor_bar(3, 3)
        real(dp) :: inverse_jacobian(3, 3), inverse_jacobian_bar(3, 3)
        real(dp) :: jacobian_inverse_bar(3, 3), determinant
        real(dp) :: determinant_bar, determinant_jacobian_bar(3, 3)
        integer :: inverse_status

        jacobian_bar = 0.0_dp
        parallel_coefficient_bar = 0.0_dp
        perpendicular_coefficient_bar = 0.0_dp
        direction_bar = 0.0_dp
        if (present(hall_coefficient_bar)) hall_coefficient_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "field-aligned tensor pullback VJP received invalid arrays")
        if (size(reference_tensor_bar, 1) /= 3 .or. &
            size(reference_tensor_bar, 2) /= 3 .or. size(direction_bar) /= 3 .or. &
            .not. all(ieee_is_finite(reference_tensor_bar)) .or. &
            .not. valid_jacobian(jacobian)) return
        call evaluate_constitutive( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            physical_tensor, status, hall_coefficient)
        if (status%code /= FORTSPARSE_OK) return
        determinant = det3(jacobian)
        call inv3(jacobian, inverse_jacobian, inverse_status)
        if (inverse_status /= 0) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "field-aligned tensor pullback VJP map is singular")
            return
        end if
        physical_tensor_bar = determinant*matmul( &
            transpose(inverse_jacobian), matmul(reference_tensor_bar, inverse_jacobian))
        inverse_jacobian_bar = determinant*( &
            matmul(reference_tensor_bar, matmul(inverse_jacobian, transpose(physical_tensor))) + &
            matmul(transpose(reference_tensor_bar), matmul(inverse_jacobian, physical_tensor)))
        determinant_bar = sum(reference_tensor_bar*matmul( &
            inverse_jacobian, matmul(physical_tensor, transpose(inverse_jacobian))))
        call inv3_vjp(jacobian, inverse_jacobian_bar, inverse_jacobian, &
            jacobian_inverse_bar, inverse_status)
        if (inverse_status /= 0) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "field-aligned tensor pullback VJP inverse failed")
            return
        end if
        call det3_vjp(jacobian, determinant_bar, determinant_jacobian_bar)
        jacobian_bar = jacobian_inverse_bar + determinant_jacobian_bar
        call evaluate_constitutive_vjp( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            physical_tensor_bar, parallel_coefficient_bar, &
            perpendicular_coefficient_bar, direction_bar, status, hall_coefficient, &
            hall_coefficient_bar)
        if (status%code /= FORTSPARSE_OK) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine pullback_field_aligned_tensor_vjp

    pure subroutine evaluate_constitutive( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            tensor, status, hall_coefficient)
        real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient
        real(dp), intent(in) :: unit_direction(:)
        real(dp), intent(out) :: tensor(3, 3)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: hall_coefficient

        if (present(hall_coefficient)) then
            call evaluate_field_aligned_constitutive_tensor( &
                parallel_coefficient, perpendicular_coefficient, unit_direction, &
                tensor, status, hall_coefficient)
        else
            call evaluate_field_aligned_constitutive_tensor( &
                parallel_coefficient, perpendicular_coefficient, unit_direction, &
                tensor, status)
        end if
    end subroutine evaluate_constitutive

    pure subroutine evaluate_constitutive_jvp( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            parallel_coefficient_dot, perpendicular_coefficient_dot, direction_dot, &
            tensor_dot, status, hall_coefficient, hall_coefficient_dot)
        real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient
        real(dp), intent(in) :: unit_direction(:), direction_dot(:)
        real(dp), intent(in) :: parallel_coefficient_dot, perpendicular_coefficient_dot
        real(dp), intent(out) :: tensor_dot(3, 3)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: hall_coefficient, hall_coefficient_dot

        if (present(hall_coefficient) .and. present(hall_coefficient_dot)) then
            call evaluate_field_aligned_constitutive_tensor_jvp( &
                parallel_coefficient, perpendicular_coefficient, unit_direction, &
                parallel_coefficient_dot, perpendicular_coefficient_dot, direction_dot, &
                tensor_dot, status, hall_coefficient, hall_coefficient_dot)
        else if (present(hall_coefficient)) then
            call evaluate_field_aligned_constitutive_tensor_jvp( &
                parallel_coefficient, perpendicular_coefficient, unit_direction, &
                parallel_coefficient_dot, perpendicular_coefficient_dot, direction_dot, &
                tensor_dot, status, hall_coefficient)
        else if (present(hall_coefficient_dot)) then
            call evaluate_field_aligned_constitutive_tensor_jvp( &
                parallel_coefficient, perpendicular_coefficient, unit_direction, &
                parallel_coefficient_dot, perpendicular_coefficient_dot, direction_dot, &
                tensor_dot, status, hall_coefficient_dot=hall_coefficient_dot)
        else
            call evaluate_field_aligned_constitutive_tensor_jvp( &
                parallel_coefficient, perpendicular_coefficient, unit_direction, &
                parallel_coefficient_dot, perpendicular_coefficient_dot, direction_dot, &
                tensor_dot, status)
        end if
    end subroutine evaluate_constitutive_jvp

    pure subroutine evaluate_constitutive_vjp( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, tensor_bar, &
            parallel_coefficient_bar, perpendicular_coefficient_bar, direction_bar, &
            status, hall_coefficient, hall_coefficient_bar)
        real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient
        real(dp), intent(in) :: unit_direction(:), tensor_bar(3, 3)
        real(dp), intent(out) :: parallel_coefficient_bar, perpendicular_coefficient_bar
        real(dp), intent(out) :: direction_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: hall_coefficient
        real(dp), intent(out), optional :: hall_coefficient_bar

        if (present(hall_coefficient) .and. present(hall_coefficient_bar)) then
            call evaluate_field_aligned_constitutive_tensor_vjp( &
                parallel_coefficient, perpendicular_coefficient, unit_direction, tensor_bar, &
                parallel_coefficient_bar, perpendicular_coefficient_bar, direction_bar, &
                status, hall_coefficient, hall_coefficient_bar)
        else if (present(hall_coefficient)) then
            call evaluate_field_aligned_constitutive_tensor_vjp( &
                parallel_coefficient, perpendicular_coefficient, unit_direction, tensor_bar, &
                parallel_coefficient_bar, perpendicular_coefficient_bar, direction_bar, &
                status, hall_coefficient)
        else if (present(hall_coefficient_bar)) then
            call evaluate_field_aligned_constitutive_tensor_vjp( &
                parallel_coefficient, perpendicular_coefficient, unit_direction, tensor_bar, &
                parallel_coefficient_bar, perpendicular_coefficient_bar, direction_bar, &
                status, hall_coefficient_bar=hall_coefficient_bar)
        else
            call evaluate_field_aligned_constitutive_tensor_vjp( &
                parallel_coefficient, perpendicular_coefficient, unit_direction, tensor_bar, &
                parallel_coefficient_bar, perpendicular_coefficient_bar, direction_bar, status)
        end if
    end subroutine evaluate_constitutive_vjp

    pure logical function valid_jacobian(jacobian) result(valid)
        real(dp), intent(in) :: jacobian(3, 3)
        real(dp) :: determinant, scale

        valid = .false.
        if (.not. all(ieee_is_finite(jacobian))) return
        determinant = det3(jacobian)
        scale = max(1.0_dp, maxval(abs(jacobian)))
        valid = determinant > 64.0_dp*epsilon(1.0_dp)*scale**3
    end function valid_jacobian

end module fortfem_field_aligned_tensor_pullback
