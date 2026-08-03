program test_generated_field_aligned_hall_products
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_field_aligned_constitutive_tensor, &
        evaluate_field_aligned_constitutive_tensor_jvp, &
        evaluate_field_aligned_constitutive_tensor_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: hall = 0.7_dp
    real(dp), parameter :: parallel = 0.0_dp
    real(dp), parameter :: perpendicular = 0.0_dp
    real(dp), parameter :: direction(3) = [0.6_dp, 0.8_dp, 0.0_dp]
    real(dp), parameter :: hall_dot = -0.3_dp
    real(dp), parameter :: direction_dot(3) = [-0.8_dp, 0.6_dp, 0.0_dp]
    real(dp), parameter :: product_bar(3) = [0.6_dp, -0.7_dp, 0.9_dp]
    real(dp), parameter :: vector(3) = [1.2_dp, -0.7_dp, 0.4_dp]
    real(dp), parameter :: epsilon = 1.0e-7_dp
    real(dp) :: product(3), product_dot(3), product_plus(3), product_minus(3)
    real(dp) :: hall_bar, direction_bar(3), lhs, rhs
    real(dp) :: tensor(3, 3), tensor_dot(3, 3), tensor_plus(3, 3)
    real(dp) :: tensor_minus(3, 3), tensor_bar(3, 3)
    real(dp) :: parallel_bar, perpendicular_bar, power
    type(fortsparse_status_t) :: status

    call evaluate_field_aligned_constitutive_tensor( &
        parallel, perpendicular, direction, tensor, status, hall)
    product = [tensor(3, 2), tensor(1, 3), tensor(2, 1)]
    call check_condition(status%code == 0 .and. &
        maxval(abs(product - hall*direction)) < 1.0e-14_dp, &
        "Hall product exposed by FEEC facade matches independent scalar products")
    power = dot_product(vector, matmul(tensor, vector))
    call check_condition(abs(power) < 1.0e-14_dp, &
        "the skew Hall tensor has zero independent instantaneous power")

    call evaluate_field_aligned_constitutive_tensor_jvp( &
        parallel, perpendicular, direction, 0.0_dp, 0.0_dp, direction_dot, &
        tensor_dot, status, hall, hall_dot)
    product_dot = [tensor_dot(3, 2), tensor_dot(1, 3), tensor_dot(2, 1)]
    call evaluate_field_aligned_constitutive_tensor( &
        parallel, perpendicular, direction + epsilon*direction_dot, &
        tensor_plus, status, hall + epsilon*hall_dot)
    call evaluate_field_aligned_constitutive_tensor( &
        parallel, perpendicular, direction - epsilon*direction_dot, &
        tensor_minus, status, hall - epsilon*hall_dot)
    product_plus = [tensor_plus(3, 2), tensor_plus(1, 3), tensor_plus(2, 1)]
    product_minus = [tensor_minus(3, 2), tensor_minus(1, 3), tensor_minus(2, 1)]
    call check_condition(maxval(abs(product_dot - &
        (product_plus - product_minus)/(2.0_dp*epsilon))) < 1.0e-9_dp, &
        "Hall-product JVP matches an independent central difference")

    tensor_bar = 0.0_dp
    tensor_bar(3, 2) = 0.5_dp*product_bar(1)
    tensor_bar(2, 3) = -0.5_dp*product_bar(1)
    tensor_bar(1, 3) = 0.5_dp*product_bar(2)
    tensor_bar(3, 1) = -0.5_dp*product_bar(2)
    tensor_bar(2, 1) = 0.5_dp*product_bar(3)
    tensor_bar(1, 2) = -0.5_dp*product_bar(3)
    call evaluate_field_aligned_constitutive_tensor_vjp( &
        parallel, perpendicular, direction, tensor_bar, parallel_bar, &
        perpendicular_bar, direction_bar, status, hall, hall_bar)
    lhs = dot_product(product_bar, product_dot)
    rhs = parallel_bar*0.0_dp + perpendicular_bar*0.0_dp + &
        hall_bar*hall_dot + dot_product(direction_bar, direction_dot)
    call check_condition(abs(lhs - rhs) < 1.0e-14_dp, &
        "Hall-product VJP satisfies the real dot-product identity")

    call check_summary("field-aligned Hall products through FEEC")
end program test_generated_field_aligned_hall_products
