program test_force_balance_generated_product
    use check, only: check_condition, check_summary
    use fortfem_feec, only: evaluate_force_balance_product, &
        evaluate_force_balance_product_jvp, evaluate_force_balance_product_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: weight = 1.7_dp
    real(dp), parameter :: test_value = -0.8_dp
    real(dp), parameter :: force_value = 2.3_dp
    real(dp), parameter :: weight_dot = -0.13_dp
    real(dp), parameter :: test_value_dot = 0.21_dp
    real(dp), parameter :: force_value_dot = -0.37_dp
    real(dp), parameter :: value_bar = -0.61_dp
    real(dp), parameter :: epsilon = 1.0e-7_dp
    real(dp) :: value, value_dot, value_plus, value_minus
    real(dp) :: weight_bar, test_value_bar, force_value_bar
    real(dp) :: lhs, rhs

    call evaluate_force_balance_product( &
        weight, test_value, force_value, value)
    call check_condition( &
        abs(value - weight*test_value*force_value) < 1.0e-14_dp, &
        "generated force product matches independent scalar oracle")

    call evaluate_force_balance_product_jvp( &
        weight, test_value, force_value, weight_dot, test_value_dot, &
        force_value_dot, value_dot)
    call evaluate_force_balance_product( &
        weight + epsilon*weight_dot, test_value + epsilon*test_value_dot, &
        force_value + epsilon*force_value_dot, value_plus)
    call evaluate_force_balance_product( &
        weight - epsilon*weight_dot, test_value - epsilon*test_value_dot, &
        force_value - epsilon*force_value_dot, value_minus)
    call check_condition( &
        abs(value_dot - (value_plus - value_minus)/(2.0_dp*epsilon)) < 2.0e-8_dp, &
        "generated force-product JVP matches independent central difference")

    call evaluate_force_balance_product_vjp( &
        weight, test_value, force_value, value_bar, weight_bar, test_value_bar, &
        force_value_bar)
    lhs = value_bar*value_dot
    rhs = weight_bar*weight_dot + test_value_bar*test_value_dot + &
        force_value_bar*force_value_dot
    call check_condition( &
        abs(lhs - rhs) < 1.0e-14_dp, &
        "generated force-product VJP satisfies the real dot-product identity")

    call check_summary("generated force-balance product")
end program test_force_balance_generated_product
