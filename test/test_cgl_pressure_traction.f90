program test_cgl_pressure_traction
    use check, only: check_condition, check_summary
    use fortfem_feec, only: evaluate_cgl_pressure_traction, &
        evaluate_cgl_pressure_traction_jvp, evaluate_cgl_pressure_traction_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: p_parallel = 3.0_dp
    real(dp), parameter :: p_perpendicular = 1.0_dp
    real(dp), parameter :: direction(3) = [0.6_dp, 0.8_dp, 0.0_dp]
    real(dp), parameter :: normal(3) = [0.5_dp, -0.2_dp, 0.7_dp]
    real(dp), parameter :: p_parallel_dot = 0.2_dp
    real(dp), parameter :: p_perpendicular_dot = -0.1_dp
    real(dp), parameter :: direction_dot(3) = [0.2_dp, -0.15_dp, 0.05_dp]
    real(dp), parameter :: normal_dot(3) = [-0.1_dp, 0.03_dp, 0.08_dp]
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp), parameter :: expected_traction(3) = [ &
        0.668_dp, 0.024_dp, 0.7_dp]
    real(dp), parameter :: traction_bar(3) = [0.4_dp, -0.3_dp, 0.8_dp]
    real(dp) :: traction(3), traction_dot(3), traction_plus(3), traction_minus(3)
    real(dp) :: p_parallel_bar, p_perpendicular_bar
    real(dp) :: direction_bar(3), normal_bar(3), lhs, rhs
    type(fortsparse_status_t) :: status

    call evaluate_cgl_pressure_traction( &
        p_parallel, p_perpendicular, direction, normal, traction, status)
    call check_condition(status%code == 0, &
        "CGL pressure traction accepts a unit field direction")
    call check_condition(maxval(abs(traction - expected_traction)) < 1.0e-14_dp, &
        "CGL pressure traction matches the independent gyrotropic oracle")

    call evaluate_cgl_pressure_traction_jvp( &
        p_parallel, p_perpendicular, direction, normal, p_parallel_dot, &
        p_perpendicular_dot, direction_dot, normal_dot, traction_dot, status)
    call evaluate_cgl_pressure_traction( &
        p_parallel + step*p_parallel_dot, p_perpendicular + step*p_perpendicular_dot, &
        direction + step*direction_dot, normal + step*normal_dot, traction_plus, &
        status)
    call evaluate_cgl_pressure_traction( &
        p_parallel - step*p_parallel_dot, p_perpendicular - step*p_perpendicular_dot, &
        direction - step*direction_dot, normal - step*normal_dot, traction_minus, &
        status)
    call check_condition(maxval(abs(traction_dot - &
        (traction_plus - traction_minus)/(2.0_dp*step))) < 3.0e-9_dp, &
        "CGL pressure traction JVP matches central differences")

    call evaluate_cgl_pressure_traction_vjp( &
        p_parallel, p_perpendicular, direction, normal, traction_bar, &
        p_parallel_bar, p_perpendicular_bar, direction_bar, normal_bar, status)
    lhs = dot_product(traction_bar, traction_dot)
    rhs = p_parallel_bar*p_parallel_dot + &
        p_perpendicular_bar*p_perpendicular_dot + &
        dot_product(direction_bar, direction_dot) + dot_product(normal_bar, normal_dot)
    call check_condition(abs(lhs - rhs) < &
        2.0e-12_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "CGL pressure traction VJP satisfies the real adjoint identity")

    call evaluate_cgl_pressure_traction( &
        p_parallel, p_perpendicular, [0.6_dp, 0.8_dp, 0.1_dp], normal, traction, &
        status)
    call check_condition(status%code /= 0, &
        "CGL pressure traction rejects a non-unit field direction")
    call check_summary("CGL pressure traction")
end program test_cgl_pressure_traction
