program test_field_aligned_tensor_pullback
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        pullback_field_aligned_tensor, pullback_field_aligned_tensor_jvp, &
        pullback_field_aligned_tensor_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: jacobian(3, 3) = reshape([ &
        1.2_dp, 0.0_dp, 0.1_dp, 0.2_dp, 0.9_dp, 0.0_dp, &
        0.1_dp, 0.2_dp, 1.1_dp], [3, 3])
    real(dp), parameter :: jacobian_dot(3, 3) = reshape([ &
        0.03_dp, -0.02_dp, 0.01_dp, -0.01_dp, 0.04_dp, 0.02_dp, &
        0.02_dp, -0.03_dp, 0.05_dp], [3, 3])
    real(dp), parameter :: direction(3) = [0.6_dp, 0.8_dp, 0.0_dp]
    real(dp), parameter :: direction_dot(3) = [0.1_dp, -0.075_dp, 0.2_dp]
    real(dp), parameter :: parallel_coefficient = 8.0_dp
    real(dp), parameter :: perpendicular_coefficient = 0.25_dp
    real(dp), parameter :: hall_coefficient = 0.35_dp
    real(dp), parameter :: parallel_dot = 0.2_dp
    real(dp), parameter :: perpendicular_dot = -0.03_dp
    real(dp), parameter :: hall_dot = -0.04_dp
    real(dp), parameter :: finite_difference_step = 1.0e-7_dp
    real(dp) :: reference_tensor(3, 3), reference_tensor_dot(3, 3)
    real(dp) :: reference_plus(3, 3), reference_minus(3, 3)
    real(dp) :: tensor(3, 3), inverse(3, 3), expected(3, 3), determinant
    real(dp) :: jacobian_bar(3, 3), direction_bar(3)
    real(dp) :: parallel_bar, perpendicular_bar, hall_bar
    real(dp) :: reference_bar(3, 3), lhs, rhs
    integer :: inverse_status
    type(fortsparse_status_t) :: status

    call pullback_field_aligned_tensor( &
        jacobian, parallel_coefficient, perpendicular_coefficient, direction, &
        reference_tensor, status, hall_coefficient)
    call check_condition(status%code == 0, &
        "field-aligned tensor pullback accepts an oriented map")

    call evaluate_independent_tensor(tensor, status)
    call independent_inverse(jacobian, inverse, determinant, inverse_status)
    expected = determinant*matmul(inverse, matmul(tensor, transpose(inverse)))
    call check_condition(inverse_status == 0 .and. maxval(abs(reference_tensor - &
        expected)) < 2.0e-13_dp, &
        "metric pullback matches the independent contravariant diffusion oracle")

    call pullback_field_aligned_tensor_jvp( &
        jacobian, parallel_coefficient, perpendicular_coefficient, direction, &
        jacobian_dot, parallel_dot, perpendicular_dot, direction_dot, &
        reference_tensor_dot, status, hall_coefficient, hall_dot)
    call pullback_field_aligned_tensor( &
        jacobian + finite_difference_step*jacobian_dot, &
        parallel_coefficient + finite_difference_step*parallel_dot, &
        perpendicular_coefficient + finite_difference_step*perpendicular_dot, &
        direction + finite_difference_step*direction_dot, reference_plus, status, &
        hall_coefficient + finite_difference_step*hall_dot)
    call pullback_field_aligned_tensor( &
        jacobian - finite_difference_step*jacobian_dot, &
        parallel_coefficient - finite_difference_step*parallel_dot, &
        perpendicular_coefficient - finite_difference_step*perpendicular_dot, &
        direction - finite_difference_step*direction_dot, reference_minus, status, &
        hall_coefficient - finite_difference_step*hall_dot)
    call check_condition(maxval(abs(reference_tensor_dot - (reference_plus - &
        reference_minus)/(2.0_dp*finite_difference_step))) < 5.0e-8_dp, &
        "metric pullback JVP matches an independent central difference")

    reference_bar = reshape([ &
        0.4_dp, -0.2_dp, 0.1_dp, 0.3_dp, 0.5_dp, -0.6_dp, &
        -0.4_dp, 0.7_dp, 0.2_dp], [3, 3])
    call pullback_field_aligned_tensor_vjp( &
        jacobian, parallel_coefficient, perpendicular_coefficient, direction, &
        reference_bar, jacobian_bar, parallel_bar, perpendicular_bar, &
        direction_bar, status, hall_coefficient, hall_bar)
    lhs = sum(reference_bar*reference_tensor_dot)
    rhs = sum(jacobian_bar*jacobian_dot) + parallel_bar*parallel_dot + &
        perpendicular_bar*perpendicular_dot + dot_product(direction_bar, direction_dot) + &
        hall_bar*hall_dot
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 3.0e-12_dp, &
        "metric pullback VJP satisfies the real transpose oracle")

    call pullback_field_aligned_tensor( &
        -jacobian, parallel_coefficient, perpendicular_coefficient, direction, &
        reference_tensor, status, hall_coefficient)
    call check_condition(status%code /= 0, &
        "metric pullback rejects an orientation-reversing map")
    call check_summary("field-aligned tensor pullback")

contains

    subroutine evaluate_independent_tensor(value, local_status)
        use fortfem_feec, only: evaluate_field_aligned_constitutive_tensor

        real(dp), intent(out) :: value(3, 3)
        type(fortsparse_status_t), intent(out) :: local_status

        call evaluate_field_aligned_constitutive_tensor( &
            parallel_coefficient, perpendicular_coefficient, direction, value, &
            local_status, hall_coefficient)
    end subroutine evaluate_independent_tensor

    subroutine independent_inverse(matrix, inverse_matrix, det, info)
        real(dp), intent(in) :: matrix(3, 3)
        real(dp), intent(out) :: inverse_matrix(3, 3), det
        integer, intent(out) :: info

        det = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - &
            matrix(2, 3)*matrix(3, 2)) - matrix(1, 2)*( &
            matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) + &
            matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - &
            matrix(2, 2)*matrix(3, 1))
        if (abs(det) < 1.0e-14_dp) then
            inverse_matrix = 0.0_dp
            info = 1
            return
        end if
        inverse_matrix(1, 1) = (matrix(2, 2)*matrix(3, 3) - &
            matrix(2, 3)*matrix(3, 2))/det
        inverse_matrix(1, 2) = (matrix(1, 3)*matrix(3, 2) - &
            matrix(1, 2)*matrix(3, 3))/det
        inverse_matrix(1, 3) = (matrix(1, 2)*matrix(2, 3) - &
            matrix(1, 3)*matrix(2, 2))/det
        inverse_matrix(2, 1) = (matrix(2, 3)*matrix(3, 1) - &
            matrix(2, 1)*matrix(3, 3))/det
        inverse_matrix(2, 2) = (matrix(1, 1)*matrix(3, 3) - &
            matrix(1, 3)*matrix(3, 1))/det
        inverse_matrix(2, 3) = (matrix(1, 3)*matrix(2, 1) - &
            matrix(1, 1)*matrix(2, 3))/det
        inverse_matrix(3, 1) = (matrix(2, 1)*matrix(3, 2) - &
            matrix(2, 2)*matrix(3, 1))/det
        inverse_matrix(3, 2) = (matrix(1, 2)*matrix(3, 1) - &
            matrix(1, 1)*matrix(3, 2))/det
        inverse_matrix(3, 3) = (matrix(1, 1)*matrix(2, 2) - &
            matrix(1, 2)*matrix(2, 1))/det
        info = 0
    end subroutine independent_inverse

end program test_field_aligned_tensor_pullback
