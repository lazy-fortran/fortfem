program test_sphere_curved_panel_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_sphere_curved_panel, evaluate_sphere_curved_panel_jvp, &
        evaluate_sphere_curved_panel_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: radius = 2.0_dp, step = 2.0e-7_dp
    real(dp) :: eta, eta_bar, eta_dot, eta_minus, eta_plus
    real(dp) :: jacobian, jacobian_bar, jacobian_dot
    real(dp) :: jacobian_minus, jacobian_plus, radius_bar, radius_dot
    real(dp) :: vertices(3, 3), vertices_bar(3, 3), vertices_dot(3, 3)
    real(dp) :: point(3), point_bar(3), point_dot(3)
    real(dp) :: point_minus(3), point_plus(3)
    real(dp) :: tangent_eta(3), tangent_eta_bar(3), tangent_eta_dot(3)
    real(dp) :: tangent_eta_minus(3), tangent_eta_plus(3)
    real(dp) :: tangent_xi(3), tangent_xi_bar(3), tangent_xi_dot(3)
    real(dp) :: tangent_xi_minus(3), tangent_xi_plus(3)
    real(dp) :: xi, xi_bar, xi_dot, xi_minus, xi_plus, lhs, rhs
    integer :: status, status_minus, status_plus

    vertices = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 0.0_dp, 1.0_dp, &
        0.0_dp, 1.0_dp, 1.0_dp], [3, 3])
    vertices_dot = reshape([ &
        0.018_dp, -0.011_dp, 0.009_dp, &
        -0.013_dp, 0.021_dp, 0.016_dp, &
        0.007_dp, -0.019_dp, 0.014_dp], [3, 3])
    xi = 0.23_dp
    eta = 0.31_dp
    xi_dot = -0.04_dp
    eta_dot = 0.025_dp
    radius_dot = 0.03_dp
    point_bar = [0.3_dp, -0.2_dp, 0.4_dp]
    tangent_xi_bar = [-0.11_dp, 0.07_dp, 0.19_dp]
    tangent_eta_bar = [0.08_dp, 0.16_dp, -0.13_dp]
    jacobian_bar = -0.35_dp

    call evaluate_sphere_curved_panel( &
        vertices, radius, xi, eta, point, tangent_xi, tangent_eta, jacobian, &
        status)
    call check_condition(status == 0, "sphere panel geometry succeeds")
    call check_condition(abs(norm2(point) - radius) < 2.0e-14_dp, &
        "sphere panel point satisfies the exact radius constraint")

    call evaluate_sphere_curved_panel_jvp( &
        vertices, radius, xi, eta, vertices_dot, radius_dot, xi_dot, eta_dot, &
        point_dot, tangent_xi_dot, tangent_eta_dot, jacobian_dot, status)
    xi_plus = xi + step*xi_dot
    eta_plus = eta + step*eta_dot
    xi_minus = xi - step*xi_dot
    eta_minus = eta - step*eta_dot
    call evaluate_sphere_curved_panel( &
        vertices + step*vertices_dot, radius + step*radius_dot, xi_plus, &
        eta_plus, point_plus, tangent_xi_plus, tangent_eta_plus, &
        jacobian_plus, status_plus)
    call evaluate_sphere_curved_panel( &
        vertices - step*vertices_dot, radius - step*radius_dot, xi_minus, &
        eta_minus, point_minus, tangent_xi_minus, tangent_eta_minus, &
        jacobian_minus, status_minus)
    call check_condition(status == 0, "sphere panel JVP succeeds")
    call check_condition(status_plus == 0, "positive sphere panel perturbation succeeds")
    call check_condition(status_minus == 0, "negative sphere panel perturbation succeeds")
    call check_condition(maxval(abs( &
        point_dot - (point_plus - point_minus)/(2.0_dp*step))) < 3.0e-8_dp, &
        "sphere point JVP matches central differences")
    call check_condition(maxval(abs( &
        tangent_xi_dot - (tangent_xi_plus - tangent_xi_minus)/(2.0_dp*step))) < &
        3.0e-8_dp, "sphere xi tangent JVP matches central differences")
    call check_condition(maxval(abs( &
        tangent_eta_dot - (tangent_eta_plus - tangent_eta_minus)/(2.0_dp*step))) < &
        3.0e-8_dp, "sphere eta tangent JVP matches central differences")
    call check_condition(abs( &
        jacobian_dot - (jacobian_plus - jacobian_minus)/(2.0_dp*step)) < 3.0e-8_dp, &
        "sphere surface Jacobian JVP matches central differences")

    call evaluate_sphere_curved_panel_vjp( &
        vertices, radius, xi, eta, point_bar, tangent_xi_bar, tangent_eta_bar, &
        jacobian_bar, vertices_bar, radius_bar, xi_bar, eta_bar, status)
    lhs = dot_product(point_bar, point_dot) + &
        dot_product(tangent_xi_bar, tangent_xi_dot) + &
        dot_product(tangent_eta_bar, tangent_eta_dot) + jacobian_bar*jacobian_dot
    rhs = sum(vertices_bar*vertices_dot) + radius_bar*radius_dot + &
        xi_bar*xi_dot + eta_bar*eta_dot
    call check_condition(status == 0, "sphere panel VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 5.0e-12_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "sphere geometry products obey the adjoint identity")

    call evaluate_sphere_curved_panel( &
        vertices, 0.0_dp, xi, eta, point, tangent_xi, tangent_eta, jacobian, &
        status)
    call check_condition(status /= 0, "sphere panel rejects a nonpositive radius")
    call evaluate_sphere_curved_panel_vjp( &
        vertices, radius, -0.1_dp, eta, point_bar, tangent_xi_bar, &
        tangent_eta_bar, jacobian_bar, vertices_bar, radius_bar, xi_bar, eta_bar, &
        status)
    call check_condition(status /= 0, "sphere panel rejects an invalid reference point")

    call check_summary("Differentiable curved sphere panel geometry")
end program test_sphere_curved_panel_ad
