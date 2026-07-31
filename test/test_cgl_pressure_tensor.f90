program test_cgl_pressure_tensor
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_cgl_pressure_tensor, evaluate_cgl_pressure_tensor_jvp, &
        evaluate_cgl_pressure_tensor_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: p_parallel = 3.0_dp
    real(dp), parameter :: p_perpendicular = 1.0_dp
    real(dp), parameter :: direction(3) = [0.6_dp, 0.8_dp, 0.0_dp]
    real(dp), parameter :: p_parallel_dot = 0.2_dp
    real(dp), parameter :: p_perpendicular_dot = -0.1_dp
    real(dp), parameter :: direction_dot(3) = [0.2_dp, -0.15_dp, 0.05_dp]
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp), parameter :: expected_tensor(3, 3) = reshape([ &
        1.72_dp, 0.96_dp, 0.0_dp, 0.96_dp, 2.28_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
    real(dp), parameter :: tensor_bar(3, 3) = reshape([ &
        0.3_dp, -0.2_dp, 0.4_dp, -0.2_dp, 0.5_dp, -0.6_dp, &
        0.4_dp, -0.6_dp, 0.7_dp], [3, 3])
    real(dp) :: pressure_tensor(3, 3), pressure_tensor_dot(3, 3)
    real(dp) :: tensor_plus(3, 3), tensor_minus(3, 3)
    real(dp) :: p_parallel_bar, p_perpendicular_bar, direction_bar(3)
    real(dp) :: vjp_left, vjp_right
    type(fortsparse_status_t) :: status

    call evaluate_cgl_pressure_tensor( &
        p_parallel, p_perpendicular, direction, pressure_tensor, status)
    call check_condition(status%code == 0, &
        "CGL pressure tensor accepts a unit magnetic direction")
    call check_condition(maxval(abs(pressure_tensor - expected_tensor)) < &
        1.0e-14_dp, &
        "CGL pressure tensor matches the independent gyrotropic oracle")

    call evaluate_cgl_pressure_tensor_jvp( &
        p_parallel, p_perpendicular, direction, p_parallel_dot, &
        p_perpendicular_dot, direction_dot, pressure_tensor_dot, status)
    call evaluate_cgl_pressure_tensor( &
        p_parallel + step*p_parallel_dot, &
        p_perpendicular + step*p_perpendicular_dot, &
        direction + step*direction_dot, tensor_plus, status)
    call evaluate_cgl_pressure_tensor( &
        p_parallel - step*p_parallel_dot, &
        p_perpendicular - step*p_perpendicular_dot, &
        direction - step*direction_dot, tensor_minus, status)
    call check_condition(maxval(abs( &
        (tensor_plus - tensor_minus)/(2.0_dp*step) - pressure_tensor_dot)) < &
        2.0e-8_dp, &
        "CGL pressure tensor JVP matches an independent central difference")

    call evaluate_cgl_pressure_tensor_vjp( &
        p_parallel, p_perpendicular, direction, tensor_bar, p_parallel_bar, &
        p_perpendicular_bar, direction_bar, status)
    vjp_left = sum(tensor_bar*pressure_tensor_dot)
    vjp_right = p_parallel_bar*p_parallel_dot + &
        p_perpendicular_bar*p_perpendicular_dot + &
        dot_product(direction_bar, direction_dot)
    call check_condition(abs(vjp_left - vjp_right) < 2.0e-13_dp, &
        "CGL pressure tensor VJP satisfies the real adjoint identity")

    call evaluate_cgl_pressure_tensor( &
        p_parallel, p_perpendicular, [1.0_dp, 0.0_dp, 1.0_dp], &
        pressure_tensor, status)
    call check_condition(status%code /= 0, &
        "CGL pressure tensor rejects a non-unit direction")
    call check_summary("Tensor-valued CGL pressure")
end program test_cgl_pressure_tensor
