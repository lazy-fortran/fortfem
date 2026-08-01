program test_interface_vector_traction_balance
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_traction_jump, &
        assemble_traction_jump_jvp, assemble_traction_jump_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: plus_traction(3) = [2.0_dp, -1.0_dp, 3.0_dp]
    real(dp), parameter :: minus_traction(3) = [0.5_dp, 0.5_dp, 1.0_dp]
    real(dp), parameter :: target(3) = [0.25_dp, -0.75_dp, 0.5_dp]
    real(dp), parameter :: plus_dot(3) = [0.2_dp, -0.3_dp, 0.4_dp]
    real(dp), parameter :: minus_dot(3) = [-0.1_dp, 0.5_dp, -0.2_dp]
    real(dp), parameter :: target_dot(3) = [0.05_dp, -0.15_dp, 0.1_dp]
    real(dp), parameter :: residual_bar(3) = [0.6_dp, -0.3_dp, 0.8_dp]
    real(dp) :: residual(3), residual_dot(3)
    real(dp) :: plus_bar(3), minus_bar(3), target_bar(3)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_traction_jump( &
        plus_traction, minus_traction, target, residual, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(residual - [1.25_dp, -0.75_dp, 1.5_dp])) < 1.0e-14_dp, &
        "vector traction jump matches the independent component oracle")

    call assemble_traction_jump_jvp( &
        plus_traction, minus_traction, target, plus_dot, minus_dot, target_dot, &
        residual_dot, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(residual_dot - [0.25_dp, -0.65_dp, 0.5_dp])) < 1.0e-14_dp, &
        "vector traction jump JVP matches the product-rule oracle")

    call assemble_traction_jump_vjp( &
        plus_traction, minus_traction, target, residual_bar, plus_bar, &
        minus_bar, target_bar, status)
    lhs = dot_product(residual_bar, residual_dot)
    rhs = dot_product(plus_bar, plus_dot) + dot_product(minus_bar, minus_dot) + &
        dot_product(target_bar, target_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-14_dp, &
        "vector traction jump VJP satisfies the real dot-product identity")

    call assemble_traction_jump( &
        plus_traction, minus_traction, target(:2), residual, status)
    call check_condition(status%code /= 0, &
        "vector traction jump rejects an incompatible target")
    call check_summary("interface vector traction balance")
end program test_interface_vector_traction_balance
