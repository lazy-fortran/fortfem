program test_cgl_pressure_work
    use check, only: check_condition, check_summary
    use fortfem_feec, only: evaluate_cgl_pressure_work, &
        evaluate_cgl_pressure_work_jvp, evaluate_cgl_pressure_work_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: p_parallel = 3.0_dp
    real(dp), parameter :: p_perpendicular = 1.0_dp
    real(dp), parameter :: direction(3) = [0.6_dp, 0.8_dp, 0.0_dp]
    real(dp), parameter :: velocity_gradient(3, 3) = reshape([ &
        0.2_dp, -0.4_dp, 0.7_dp, 0.1_dp, 0.3_dp, -0.2_dp, &
        0.5_dp, 0.6_dp, -0.1_dp], [3, 3])
    real(dp), parameter :: p_parallel_dot = 0.2_dp
    real(dp), parameter :: p_perpendicular_dot = -0.1_dp
    real(dp), parameter :: direction_dot(3) = [0.08_dp, -0.06_dp, 0.03_dp]
    real(dp), parameter :: velocity_gradient_dot(3, 3) = reshape([ &
        -0.3_dp, 0.2_dp, 0.1_dp, 0.4_dp, -0.2_dp, 0.5_dp, &
        -0.1_dp, 0.3_dp, 0.6_dp], [3, 3])
    real(dp), parameter :: work_bar = 0.7_dp
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: work, work_dot, work_plus, work_minus
    real(dp) :: p_parallel_bar, p_perpendicular_bar, direction_bar(3)
    real(dp) :: velocity_gradient_bar(3, 3), lhs, rhs
    real(dp) :: direction_plus(3), direction_minus(3)
    real(dp) :: bad_gradient(2, 2)
    type(fortsparse_status_t) :: status

    call evaluate_cgl_pressure_work( &
        p_parallel, p_perpendicular, direction, velocity_gradient, work, status)
    call check_condition(status%code == 0, &
        "CGL pressure work accepts a unit field direction")
    call check_condition(abs(work - 0.64_dp) < 1.0e-14_dp, &
        "CGL pressure work matches the independent tensor contraction oracle")

    call evaluate_cgl_pressure_work_jvp( &
        p_parallel, p_perpendicular, direction, velocity_gradient, &
        p_parallel_dot, p_perpendicular_dot, direction_dot, &
        velocity_gradient_dot, work_dot, status)
    direction_plus = direction + step*direction_dot
    direction_minus = direction - step*direction_dot
    call evaluate_cgl_pressure_work( &
        p_parallel + step*p_parallel_dot, &
        p_perpendicular + step*p_perpendicular_dot, direction_plus, &
        velocity_gradient + step*velocity_gradient_dot, work_plus, status)
    call evaluate_cgl_pressure_work( &
        p_parallel - step*p_parallel_dot, &
        p_perpendicular - step*p_perpendicular_dot, direction_minus, &
        velocity_gradient - step*velocity_gradient_dot, work_minus, status)
    call check_condition(abs(work_dot - (work_plus - work_minus)/(2.0_dp*step)) &
        < 3.0e-9_dp, "CGL pressure work JVP matches central differences")

    call evaluate_cgl_pressure_work_vjp( &
        p_parallel, p_perpendicular, direction, velocity_gradient, work_bar, &
        p_parallel_bar, p_perpendicular_bar, direction_bar, &
        velocity_gradient_bar, status)
    lhs = work_bar*work_dot
    rhs = p_parallel_bar*p_parallel_dot + &
        p_perpendicular_bar*p_perpendicular_dot + &
        dot_product(direction_bar, direction_dot) + &
        sum(velocity_gradient_bar*velocity_gradient_dot)
    call check_condition(abs(lhs - rhs) < &
        2.0e-12_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "CGL pressure work VJP satisfies the real adjoint identity")

    bad_gradient = 0.0_dp
    call evaluate_cgl_pressure_work( &
        p_parallel, p_perpendicular, direction, bad_gradient, work, status)
    call check_condition(status%code /= 0, &
        "CGL pressure work rejects an incompatible gradient")
    call check_summary("CGL pressure work")
end program test_cgl_pressure_work
