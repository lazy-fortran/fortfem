program test_surface_triangle_geometry_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_surface_triangle_geometry_3d, &
        evaluate_surface_triangle_geometry_3d_jvp, &
        evaluate_surface_triangle_geometry_3d_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: eta, eta_bar, eta_dot, eta_minus, eta_plus
    real(dp) :: jacobian, jacobian_bar, jacobian_dot
    real(dp) :: jacobian_minus, jacobian_plus
    real(dp) :: normal(3), normal_bar(3), normal_dot(3)
    real(dp) :: normal_minus(3), normal_plus(3)
    real(dp) :: point(3), point_bar(3), point_dot(3)
    real(dp) :: point_minus(3), point_plus(3)
    real(dp) :: vertices(3, 3), vertices_bar(3, 3), vertices_dot(3, 3)
    real(dp) :: xi, xi_bar, xi_dot, xi_minus, xi_plus
    real(dp) :: lhs, rhs
    integer :: status, status_minus, status_plus

    vertices = reshape([ &
        0.2_dp, -0.1_dp, 0.3_dp, &
        1.4_dp, 0.2_dp, -0.15_dp, &
        -0.25_dp, 1.1_dp, 0.6_dp], [3, 3])
    vertices_dot = reshape([ &
        0.03_dp, -0.02_dp, 0.01_dp, &
        -0.015_dp, 0.04_dp, 0.025_dp, &
        0.02_dp, -0.01_dp, 0.035_dp], [3, 3])
    xi = 0.23_dp
    eta = 0.31_dp
    xi_dot = -0.04_dp
    eta_dot = 0.025_dp
    point_bar = [0.3_dp, -0.2_dp, 0.4_dp]
    jacobian_bar = -0.35_dp
    normal_bar = [-0.1_dp, 0.25_dp, 0.15_dp]

    call evaluate_surface_triangle_geometry_3d( &
        vertices, xi, eta, point, jacobian, normal, status)
    call check_condition(status == 0, "surface triangle geometry succeeds")
    call check_condition(abs(jacobian - 1.703622977656735_dp) < 2.0e-14_dp, &
        "surface Jacobian matches an independent numerical oracle")
    call check_condition(maxval(abs( &
        point - [0.3365_dp, 0.341_dp, 0.2895_dp])) < 2.0e-15_dp, &
        "affine surface point matches barycentric evaluation")

    call evaluate_surface_triangle_geometry_3d_jvp( &
        vertices, xi, eta, vertices_dot, xi_dot, eta_dot, point_dot, &
        jacobian_dot, normal_dot, status)
    xi_plus = xi + step*xi_dot
    eta_plus = eta + step*eta_dot
    xi_minus = xi - step*xi_dot
    eta_minus = eta - step*eta_dot
    call evaluate_surface_triangle_geometry_3d( &
        vertices + step*vertices_dot, xi_plus, eta_plus, point_plus, &
        jacobian_plus, normal_plus, status_plus)
    call evaluate_surface_triangle_geometry_3d( &
        vertices - step*vertices_dot, xi_minus, eta_minus, point_minus, &
        jacobian_minus, normal_minus, status_minus)
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "surface geometry JVP accepts a valid direction")
    call check_condition(maxval(abs( &
        point_dot - (point_plus - point_minus)/(2.0_dp*step))) < 2.0e-9_dp, &
        "surface point JVP matches central differences")
    call check_condition(abs( &
        jacobian_dot - (jacobian_plus - jacobian_minus)/(2.0_dp*step)) < &
        2.0e-9_dp, "surface Jacobian JVP matches central differences")
    call check_condition(maxval(abs( &
        normal_dot - (normal_plus - normal_minus)/(2.0_dp*step))) < &
        2.0e-9_dp, "surface normal JVP matches central differences")

    call evaluate_surface_triangle_geometry_3d_vjp( &
        vertices, xi, eta, point_bar, jacobian_bar, normal_bar, vertices_bar, &
        xi_bar, eta_bar, status)
    lhs = dot_product(point_bar, point_dot) + jacobian_bar*jacobian_dot + &
        dot_product(normal_bar, normal_dot)
    rhs = sum(vertices_bar*vertices_dot) + xi_bar*xi_dot + eta_bar*eta_dot
    call check_condition(status == 0, "surface geometry VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 3.0e-13_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "surface geometry products obey the adjoint identity")

    vertices(:, 3) = vertices(:, 2)
    call evaluate_surface_triangle_geometry_3d( &
        vertices, xi, eta, point, jacobian, normal, status)
    call check_condition( &
        status /= 0 .and. jacobian == 0.0_dp .and. all(normal == 0.0_dp), &
        "surface geometry rejects a degenerate triangle")

    call check_summary("Surface triangle geometry derivatives")
end program test_surface_triangle_geometry_ad
