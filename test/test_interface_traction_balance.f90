program test_interface_traction_balance
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_normal_traction_jump, &
        assemble_normal_traction_jump_jvp, assemble_normal_traction_jump_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: plus_traction(3) = [2.0_dp, -1.0_dp, 3.0_dp]
    real(dp), parameter :: minus_traction(3) = [0.5_dp, 0.5_dp, 1.0_dp]
    real(dp), parameter :: normal(3) = [0.0_dp, 0.0_dp, 1.0_dp]
    real(dp), parameter :: target = 0.75_dp
    real(dp), parameter :: plus_dot(3) = [0.2_dp, -0.3_dp, 0.4_dp]
    real(dp), parameter :: minus_dot(3) = [-0.1_dp, 0.5_dp, -0.2_dp]
    real(dp), parameter :: normal_dot(3) = [0.1_dp, -0.2_dp, 0.0_dp]
    real(dp), parameter :: target_dot = -0.15_dp
    real(dp), parameter :: residual_bar = 0.6_dp
    real(dp) :: residual, residual_dot
    real(dp) :: plus_bar(3), minus_bar(3), normal_bar(3), target_bar
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_normal_traction_jump( &
        plus_traction, minus_traction, normal, target, residual, status)
    call check_condition(status%code == 0 .and. &
        abs(residual - 1.25_dp) < 1.0e-14_dp, &
        "normal traction jump matches the independent scalar oracle")

    call assemble_normal_traction_jump_jvp( &
        plus_traction, minus_traction, normal, target, plus_dot, minus_dot, &
        normal_dot, target_dot, residual_dot, status)
    call check_condition(status%code == 0 .and. &
        abs(residual_dot - 1.2_dp) < 1.0e-14_dp, &
        "normal traction jump JVP matches the product-rule oracle")

    call assemble_normal_traction_jump_vjp( &
        plus_traction, minus_traction, normal, target, residual_bar, plus_bar, &
        minus_bar, normal_bar, target_bar, status)
    lhs = residual_bar*residual_dot
    rhs = dot_product(plus_bar, plus_dot) + dot_product(minus_bar, minus_dot) + &
        dot_product(normal_bar, normal_dot) + target_bar*target_dot
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-14_dp, &
        "normal traction jump VJP satisfies the real dot-product identity")

    call assemble_normal_traction_jump( &
        plus_traction, minus_traction, [0.0_dp, 0.0_dp, 2.0_dp], target, &
        residual, status)
    call check_condition(status%code /= 0, &
        "normal traction jump rejects a non-unit normal")

    call check_summary("interface normal traction balance")
end program test_interface_traction_balance
