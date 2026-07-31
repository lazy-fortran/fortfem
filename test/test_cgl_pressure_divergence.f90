program test_cgl_pressure_divergence
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_cgl_pressure_divergence, &
        evaluate_cgl_pressure_divergence_jvp, &
        evaluate_cgl_pressure_divergence_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: p_parallel = 3.0_dp
    real(dp), parameter :: p_perpendicular = 1.0_dp
    real(dp), parameter :: direction(3) = [1.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: parallel_gradient(3) = [2.0_dp, 3.0_dp, 4.0_dp]
    real(dp), parameter :: perpendicular_gradient(3) = [1.0_dp, -1.0_dp, 0.0_dp]
    real(dp), parameter :: direction_gradient(3, 3) = reshape([ &
        0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp], [3, 3])
    real(dp), parameter :: expected_force(3) = [2.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: p_parallel_dot = 0.2_dp
    real(dp), parameter :: p_perpendicular_dot = -0.1_dp
    real(dp), parameter :: direction_dot(3) = [0.0_dp, 0.1_dp, 0.0_dp]
    real(dp), parameter :: parallel_gradient_dot(3) = [0.1_dp, -0.2_dp, 0.3_dp]
    real(dp), parameter :: perpendicular_gradient_dot(3) = [-0.2_dp, 0.3_dp, 0.1_dp]
    real(dp), parameter :: direction_gradient_dot(3, 3) = reshape([ &
        0.04_dp, -0.03_dp, 0.02_dp, 0.01_dp, 0.05_dp, -0.04_dp, &
        -0.02_dp, 0.03_dp, 0.06_dp], [3, 3])
    real(dp), parameter :: force_bar(3) = [0.3_dp, -0.2_dp, 0.4_dp]
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: force(3), force_dot(3), force_plus(3), force_minus(3)
    real(dp) :: p_parallel_bar, p_perpendicular_bar, direction_bar(3)
    real(dp) :: parallel_gradient_bar(3), perpendicular_gradient_bar(3)
    real(dp) :: direction_gradient_bar(3, 3)
    real(dp) :: vjp_left, vjp_right
    type(fortsparse_status_t) :: status

    call evaluate_cgl_pressure_divergence( &
        p_parallel, p_perpendicular, direction, parallel_gradient, &
        perpendicular_gradient, direction_gradient, force, status)
    call check_condition(status%code == 0, &
        "CGL pressure divergence accepts a unit magnetic direction")
    call check_condition(maxval(abs(force - expected_force)) < 1.0e-14_dp, &
        "CGL pressure divergence matches the independent product-rule oracle")

    call evaluate_cgl_pressure_divergence_jvp( &
        p_parallel, p_perpendicular, direction, parallel_gradient, &
        perpendicular_gradient, direction_gradient, p_parallel_dot, &
        p_perpendicular_dot, direction_dot, parallel_gradient_dot, &
        perpendicular_gradient_dot, direction_gradient_dot, force_dot, status)
    call evaluate_cgl_pressure_divergence( &
        p_parallel + step*p_parallel_dot, p_perpendicular + step*p_perpendicular_dot, &
        direction + step*direction_dot, parallel_gradient + step*parallel_gradient_dot, &
        perpendicular_gradient + step*perpendicular_gradient_dot, &
        direction_gradient + step*direction_gradient_dot, force_plus, status)
    call evaluate_cgl_pressure_divergence( &
        p_parallel - step*p_parallel_dot, p_perpendicular - step*p_perpendicular_dot, &
        direction - step*direction_dot, parallel_gradient - step*parallel_gradient_dot, &
        perpendicular_gradient - step*perpendicular_gradient_dot, &
        direction_gradient - step*direction_gradient_dot, force_minus, status)
    call check_condition(maxval(abs( &
        (force_plus - force_minus)/(2.0_dp*step) - force_dot)) < 2.0e-7_dp, &
        "CGL pressure divergence JVP matches an independent central difference")

    call evaluate_cgl_pressure_divergence_vjp( &
        p_parallel, p_perpendicular, direction, parallel_gradient, &
        perpendicular_gradient, direction_gradient, force_bar, p_parallel_bar, &
        p_perpendicular_bar, direction_bar, parallel_gradient_bar, &
        perpendicular_gradient_bar, direction_gradient_bar, status)
    vjp_left = dot_product(force_bar, force_dot)
    vjp_right = p_parallel_bar*p_parallel_dot + &
        p_perpendicular_bar*p_perpendicular_dot + dot_product(direction_bar, direction_dot) + &
        dot_product(parallel_gradient_bar, parallel_gradient_dot) + &
        dot_product(perpendicular_gradient_bar, perpendicular_gradient_dot) + &
        sum(direction_gradient_bar*direction_gradient_dot)
    call check_condition(abs(vjp_left - vjp_right) < 2.0e-12_dp, &
        "CGL pressure divergence VJP satisfies the real adjoint identity")

    call evaluate_cgl_pressure_divergence( &
        p_parallel, p_perpendicular, [1.0_dp, 0.0_dp, 1.0_dp], &
        parallel_gradient, perpendicular_gradient, direction_gradient, force, status)
    call check_condition(status%code /= 0, &
        "CGL pressure divergence rejects a non-unit direction")
    call check_summary("CGL pressure divergence")
end program test_cgl_pressure_divergence
