program test_generated_field_aligned_hall_products
    use check, only: check_condition, check_summary
    use fortfem_generated_field_aligned_hall, only: &
        generated_field_aligned_hall
    use fortfem_generated_field_aligned_hall_jvp, only: &
        generated_field_aligned_hall_jvp
    use fortfem_generated_field_aligned_hall_vjp, only: &
        generated_field_aligned_hall_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: hall = 0.7_dp
    real(dp), parameter :: direction(3) = [0.2_dp, -0.4_dp, 0.8_dp]
    real(dp), parameter :: hall_dot = -0.3_dp
    real(dp), parameter :: direction_dot(3) = [0.5_dp, 0.1_dp, -0.2_dp]
    real(dp), parameter :: product_bar(3) = [0.6_dp, -0.7_dp, 0.9_dp]
    real(dp), parameter :: epsilon = 1.0e-7_dp
    real(dp) :: product(3), product_dot(3), product_plus(3), product_minus(3)
    real(dp) :: hall_bar, direction_bar(3), lhs, rhs

    call generated_field_aligned_hall( &
        hall, direction(1), direction(2), direction(3), product(1), &
        product(2), product(3))
    call check_condition(maxval(abs(product - hall*direction)) < 1.0e-14_dp, &
        "generated Hall product matches independent scalar products")

    call generated_field_aligned_hall_jvp( &
        hall, direction(1), direction(2), direction(3), hall_dot, &
        direction_dot(1), direction_dot(2), direction_dot(3), product_dot(1), &
        product_dot(2), product_dot(3))
    call generated_field_aligned_hall( &
        hall + epsilon*hall_dot, direction(1) + epsilon*direction_dot(1), &
        direction(2) + epsilon*direction_dot(2), &
        direction(3) + epsilon*direction_dot(3), product_plus(1), &
        product_plus(2), product_plus(3))
    call generated_field_aligned_hall( &
        hall - epsilon*hall_dot, direction(1) - epsilon*direction_dot(1), &
        direction(2) - epsilon*direction_dot(2), &
        direction(3) - epsilon*direction_dot(3), product_minus(1), &
        product_minus(2), product_minus(3))
    call check_condition(maxval(abs(product_dot - &
        (product_plus - product_minus)/(2.0_dp*epsilon))) < 1.0e-9_dp, &
        "generated Hall-product JVP matches a central difference")

    call generated_field_aligned_hall_vjp( &
        hall, direction(1), direction(2), direction(3), product_bar(1), &
        product_bar(2), product_bar(3), hall_bar, direction_bar(1), &
        direction_bar(2), direction_bar(3))
    lhs = dot_product(product_bar, product_dot)
    rhs = hall_bar*hall_dot + dot_product(direction_bar, direction_dot)
    call check_condition(abs(lhs - rhs) < 1.0e-14_dp, &
        "generated Hall-product VJP satisfies the real dot-product identity")

    call check_summary("generated field-aligned Hall products")
end program test_generated_field_aligned_hall_products
