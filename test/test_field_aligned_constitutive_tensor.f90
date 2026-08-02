program test_field_aligned_constitutive_tensor
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_field_aligned_constitutive_tensor, &
        evaluate_field_aligned_constitutive_tensor_jvp, &
        evaluate_field_aligned_constitutive_tensor_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: parallel_coefficient = 3.0_dp
    real(dp), parameter :: perpendicular_coefficient = 1.0_dp
    real(dp), parameter :: hall_coefficient = 0.5_dp
    real(dp), parameter :: direction(3) = [0.6_dp, 0.8_dp, 0.0_dp]
    real(dp), parameter :: parallel_dot = 0.2_dp
    real(dp), parameter :: perpendicular_dot = -0.1_dp
    real(dp), parameter :: hall_dot = -0.3_dp
    real(dp), parameter :: direction_dot(3) = [0.2_dp, -0.15_dp, 0.1_dp]
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp), parameter :: expected_tensor(3, 3) = reshape([ &
        1.72_dp, 0.96_dp, -0.4_dp, 0.96_dp, 2.28_dp, 0.3_dp, &
        0.4_dp, -0.3_dp, 1.0_dp], [3, 3])
    real(dp), parameter :: tensor_bar(3, 3) = reshape([ &
        0.3_dp, -0.2_dp, 0.4_dp, -0.2_dp, 0.5_dp, -0.6_dp, &
        0.4_dp, -0.6_dp, 0.7_dp], [3, 3])
    real(dp) :: tensor(3, 3), tensor_dot(3, 3), tensor_plus(3, 3)
    real(dp) :: tensor_minus(3, 3), no_hall_tensor(3, 3)
    real(dp) :: parallel_bar, perpendicular_bar, hall_bar, direction_bar(3)
    real(dp) :: left, right
    type(fortsparse_status_t) :: status

    call evaluate_field_aligned_constitutive_tensor( &
        parallel_coefficient, perpendicular_coefficient, direction, tensor, &
        status, hall_coefficient)
    call check_condition(status%code == 0, &
        "field-aligned constitutive tensor accepts finite unit data")
    call check_condition(maxval(abs(tensor - expected_tensor)) < 1.0e-14_dp, &
        "symmetric projectors plus b-cross product match independent oracle")

    call evaluate_field_aligned_constitutive_tensor( &
        parallel_coefficient, perpendicular_coefficient, direction, &
        no_hall_tensor, status)
    call check_condition(status%code == 0 .and. maxval(abs(no_hall_tensor - &
        transpose(no_hall_tensor))) < 1.0e-14_dp, &
        "omitting the optional Hall coefficient gives a symmetric tensor")

    call evaluate_field_aligned_constitutive_tensor_jvp( &
        parallel_coefficient, perpendicular_coefficient, direction, parallel_dot, &
        perpendicular_dot, direction_dot, tensor_dot, status, hall_coefficient, &
        hall_dot)
    call evaluate_field_aligned_constitutive_tensor( &
        parallel_coefficient + step*parallel_dot, &
        perpendicular_coefficient + step*perpendicular_dot, &
        direction + step*direction_dot, tensor_plus, status, &
        hall_coefficient + step*hall_dot)
    call evaluate_field_aligned_constitutive_tensor( &
        parallel_coefficient - step*parallel_dot, &
        perpendicular_coefficient - step*perpendicular_dot, &
        direction - step*direction_dot, tensor_minus, status, &
        hall_coefficient - step*hall_dot)
    call check_condition(maxval(abs(tensor_dot - &
        (tensor_plus - tensor_minus)/(2.0_dp*step))) < 3.0e-8_dp, &
        "constitutive tensor JVP matches an independent central difference")

    call evaluate_field_aligned_constitutive_tensor_vjp( &
        parallel_coefficient, perpendicular_coefficient, direction, tensor_bar, &
        parallel_bar, perpendicular_bar, direction_bar, status, hall_coefficient, &
        hall_bar)
    left = sum(tensor_bar*tensor_dot)
    right = parallel_bar*parallel_dot + perpendicular_bar*perpendicular_dot + &
        hall_bar*hall_dot + dot_product(direction_bar, direction_dot)
    call check_condition(status%code == 0 .and. abs(left - right) < 2.0e-13_dp, &
        "constitutive tensor VJP satisfies the real transpose oracle")

    call evaluate_field_aligned_constitutive_tensor( &
        ieee_value(0.0_dp, ieee_quiet_nan), perpendicular_coefficient, &
        direction, tensor, status, hall_coefficient)
    call check_condition(status%code /= 0, &
        "constitutive tensor rejects non-finite coefficients")
    call evaluate_field_aligned_constitutive_tensor( &
        parallel_coefficient, perpendicular_coefficient, [0.6_dp, 0.8_dp, 0.1_dp], &
        tensor, status, hall_coefficient)
    call check_condition(status%code /= 0, &
        "constitutive tensor rejects a non-unit direction")
    call check_summary("field-aligned constitutive tensor")
end program test_field_aligned_constitutive_tensor
