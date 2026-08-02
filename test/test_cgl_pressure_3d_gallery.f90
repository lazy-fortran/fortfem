program test_cgl_pressure_3d_gallery
    !! Independent oracle for the 3-D tensor-pressure gallery contract.
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_cgl_pressure_tensor, evaluate_cgl_pressure_tensor_jvp, &
        evaluate_cgl_pressure_tensor_vjp, evaluate_cgl_pressure_traction, &
        evaluate_cgl_pressure_traction_jvp, evaluate_cgl_pressure_traction_vjp, &
        evaluate_cgl_pressure_work, evaluate_cgl_pressure_work_jvp, &
        evaluate_cgl_pressure_work_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: p_parallel = 3.7_dp
    real(dp), parameter :: p_perpendicular = 1.15_dp
    real(dp), parameter :: direction(3) = [1.0_dp, 2.0_dp, 2.0_dp]/3.0_dp
    real(dp), parameter :: direction_dot(3) = [0.04_dp, -0.03_dp, 0.02_dp]
    real(dp), parameter :: normal(3) = [0.3_dp, -0.7_dp, 0.6_dp]
    real(dp), parameter :: normal_dot(3) = [-0.02_dp, 0.05_dp, 0.03_dp]
    real(dp), parameter :: p_parallel_dot = 0.11_dp
    real(dp), parameter :: p_perpendicular_dot = -0.07_dp
    real(dp), parameter :: gradient(3, 3) = reshape([ &
        0.2_dp, -0.1_dp, 0.3_dp, 0.4_dp, 0.5_dp, -0.2_dp, &
        0.1_dp, 0.25_dp, -0.35_dp], [3, 3])
    real(dp), parameter :: gradient_dot(3, 3) = reshape([ &
        -0.3_dp, 0.15_dp, 0.2_dp, -0.1_dp, 0.25_dp, 0.4_dp, &
        0.05_dp, -0.2_dp, 0.1_dp], [3, 3])
    real(dp), parameter :: tensor_bar(3, 3) = reshape([ &
        0.2_dp, -0.4_dp, 0.5_dp, 0.1_dp, 0.7_dp, -0.3_dp, &
        -0.6_dp, 0.8_dp, 0.9_dp], [3, 3])
    real(dp), parameter :: traction_bar(3) = [0.6_dp, -0.2_dp, 0.8_dp]
    real(dp), parameter :: work_bar = 0.73_dp
    real(dp) :: tensor(3, 3), tensor_dot(3, 3), tensor_expected(3, 3)
    real(dp) :: traction(3), traction_dot(3), traction_expected(3)
    real(dp) :: work, work_dot, work_expected
    real(dp) :: p_parallel_bar, p_perpendicular_bar
    real(dp) :: direction_bar(3), normal_bar(3), gradient_bar(3, 3)
    real(dp) :: lhs, rhs, direction_dot_tangent(3)
    real(dp) :: identity(3, 3), delta, step, tensor_plus(3, 3)
    real(dp) :: tensor_minus(3, 3), traction_plus(3), traction_minus(3)
    real(dp) :: work_plus, work_minus
    type(fortsparse_status_t) :: status

    identity = 0.0_dp
    identity(1, 1) = 1.0_dp
    identity(2, 2) = 1.0_dp
    identity(3, 3) = 1.0_dp
    delta = p_parallel - p_perpendicular
    tensor_expected = p_perpendicular*identity + &
        delta*outer_product(direction, direction)

    call evaluate_cgl_pressure_tensor( &
        p_parallel, p_perpendicular, direction, tensor, status)
    call check_condition(status%code == 0, &
        "3-D CGL tensor accepts an oblique unit direction")
    call check_condition(maxval(abs(tensor - tensor_expected)) < 2.0e-14_dp, &
        "3-D tensor contraction matches the independent projector oracle")

    direction_dot_tangent = direction_dot - &
        direction*dot_product(direction, direction_dot)
    tensor_expected = p_perpendicular_dot*identity + &
        (p_parallel_dot - p_perpendicular_dot)* &
        outer_product(direction, direction) + delta*( &
        outer_product(direction_dot_tangent, direction) + &
        outer_product(direction, direction_dot_tangent))
    call evaluate_cgl_pressure_tensor_jvp( &
        p_parallel, p_perpendicular, direction, p_parallel_dot, &
        p_perpendicular_dot, direction_dot_tangent, tensor_dot, status)
    call check_condition(maxval(abs(tensor_dot - tensor_expected)) < 2.0e-13_dp, &
        "3-D tensor JVP matches the independent product-rule oracle")

    call evaluate_cgl_pressure_tensor_vjp( &
        p_parallel, p_perpendicular, direction, tensor_bar, p_parallel_bar, &
        p_perpendicular_bar, direction_bar, status)
    lhs = sum(tensor_bar*tensor_dot)
    rhs = p_parallel_bar*p_parallel_dot + &
        p_perpendicular_bar*p_perpendicular_dot + &
        dot_product(direction_bar, direction_dot_tangent)
    call check_condition(abs(lhs - rhs) < 2.0e-12_dp, &
        "3-D tensor VJP satisfies the independent adjoint identity")

    call evaluate_cgl_pressure_traction( &
        p_parallel, p_perpendicular, direction, normal, traction, status)
    traction_expected = matmul(tensor_expected + &
        (p_parallel_dot*0.0_dp)*identity, normal)
    call check_condition(status%code == 0, &
        "3-D CGL traction accepts an oblique unit direction")
    call check_condition(maxval(abs(traction - matmul(tensor, normal))) < &
        2.0e-14_dp, "traction equals the independent tensor-normal contraction")

    call evaluate_cgl_pressure_traction_jvp( &
        p_parallel, p_perpendicular, direction, normal, p_parallel_dot, &
        p_perpendicular_dot, direction_dot_tangent, normal_dot, traction_dot, &
        status)
    traction_expected = matmul(tensor_dot, normal) + matmul(tensor, normal_dot)
    call check_condition(maxval(abs(traction_dot - traction_expected)) < &
        2.0e-13_dp, "traction JVP matches the independent product rule")

    call evaluate_cgl_pressure_traction_vjp( &
        p_parallel, p_perpendicular, direction, normal, traction_bar, &
        p_parallel_bar, p_perpendicular_bar, direction_bar, normal_bar, status)
    lhs = dot_product(traction_bar, traction_dot)
    rhs = p_parallel_bar*p_parallel_dot + &
        p_perpendicular_bar*p_perpendicular_dot + &
        dot_product(direction_bar, direction_dot_tangent) + &
        dot_product(normal_bar, normal_dot)
    call check_condition(abs(lhs - rhs) < 2.0e-12_dp, &
        "traction VJP satisfies the independent adjoint identity")

    call evaluate_cgl_pressure_work( &
        p_parallel, p_perpendicular, direction, gradient, work, status)
    work_expected = sum(tensor*gradient)
    call check_condition(status%code == 0 .and. abs(work - work_expected) < &
        2.0e-14_dp, "pressure work is the independent tensor contraction")
    call evaluate_cgl_pressure_work_jvp( &
        p_parallel, p_perpendicular, direction, gradient, p_parallel_dot, &
        p_perpendicular_dot, direction_dot_tangent, gradient_dot, work_dot, &
        status)
    work_expected = sum(tensor_dot*gradient) + sum(tensor*gradient_dot)
    call check_condition(abs(work_dot - work_expected) < 2.0e-13_dp, &
        "pressure-work JVP matches the independent contraction rule")
    call evaluate_cgl_pressure_work_vjp( &
        p_parallel, p_perpendicular, direction, gradient, work_bar, &
        p_parallel_bar, p_perpendicular_bar, direction_bar, gradient_bar, &
        status)
    lhs = work_bar*work_dot
    rhs = p_parallel_bar*p_parallel_dot + &
        p_perpendicular_bar*p_perpendicular_dot + &
        dot_product(direction_bar, direction_dot_tangent) + &
        sum(gradient_bar*gradient_dot)
    call check_condition(abs(lhs - rhs) < 2.0e-12_dp, &
        "pressure-work VJP satisfies the independent adjoint identity")

    step = 1.0e-7_dp
    call evaluate_cgl_pressure_tensor( &
        p_parallel + step*p_parallel_dot, p_perpendicular + &
        step*p_perpendicular_dot, direction + step*direction_dot_tangent, &
        tensor_plus, status)
    call evaluate_cgl_pressure_tensor( &
        p_parallel - step*p_parallel_dot, p_perpendicular - &
        step*p_perpendicular_dot, direction - step*direction_dot_tangent, &
        tensor_minus, status)
    call check_condition(maxval(abs(tensor_dot - &
        (tensor_plus - tensor_minus)/(2.0_dp*step))) < 3.0e-8_dp, &
        "3-D tensor JVP agrees with an independent finite difference")

    call evaluate_cgl_pressure_traction( &
        p_parallel + step*p_parallel_dot, p_perpendicular + &
        step*p_perpendicular_dot, direction + step*direction_dot_tangent, &
        normal + step*normal_dot, traction_plus, status)
    call evaluate_cgl_pressure_traction( &
        p_parallel - step*p_parallel_dot, p_perpendicular - &
        step*p_perpendicular_dot, direction - step*direction_dot_tangent, &
        normal - step*normal_dot, traction_minus, status)
    call check_condition(maxval(abs(traction_dot - &
        (traction_plus - traction_minus)/(2.0_dp*step))) < 3.0e-8_dp, &
        "3-D traction JVP agrees with an independent finite difference")

    call evaluate_cgl_pressure_work( &
        p_parallel + step*p_parallel_dot, p_perpendicular + &
        step*p_perpendicular_dot, direction + step*direction_dot_tangent, &
        gradient + step*gradient_dot, work_plus, status)
    call evaluate_cgl_pressure_work( &
        p_parallel - step*p_parallel_dot, p_perpendicular - &
        step*p_perpendicular_dot, direction - step*direction_dot_tangent, &
        gradient - step*gradient_dot, work_minus, status)
    call check_condition(abs(work_dot - (work_plus - work_minus)/ &
        (2.0_dp*step)) < 3.0e-8_dp, &
        "pressure-work JVP agrees with an independent finite difference")
    call check_summary("3-D CGL tensor pressure gallery")

contains

    pure function outer_product(left, right) result(product)
        real(dp), intent(in) :: left(3), right(3)
        real(dp) :: product(3, 3)

        product = spread(left, dim=2, ncopies=3)* &
            spread(right, dim=1, ncopies=3)
    end function outer_product

end program test_cgl_pressure_3d_gallery
