program test_maxwell_torus_curved_rwg_basis_ad
    use check, only: check_condition, check_summary
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_feec, only: &
        build_maxwell_rwg_surface_space, &
        evaluate_maxwell_torus_curved_rwg_basis, &
        evaluate_maxwell_torus_curved_rwg_basis_jvp, &
        evaluate_maxwell_torus_curved_rwg_basis_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.1_dp, minor_radius = 0.55_dp
    real(dp), parameter :: xi = 0.23_dp, eta = 0.29_dp, step = 2.0e-6_dp
    real(dp) :: major_radius_dot, minor_radius_dot, xi_dot, eta_dot
    real(dp) :: major_radius_bar, minor_radius_bar, xi_bar, eta_bar
    real(dp) :: divergence, divergence_bar, divergence_dot
    real(dp) :: divergence_minus, divergence_plus
    real(dp) :: jacobian, jacobian_bar, jacobian_dot
    real(dp) :: jacobian_minus, jacobian_plus
    real(dp) :: lhs, rhs
    real(dp) :: point(3), point_bar(3), point_dot(3)
    real(dp) :: point_minus(3), point_plus(3)
    real(dp) :: value(3), value_bar(3), value_dot(3)
    real(dp) :: value_minus(3), value_plus(3)
    real(dp), allocatable :: parameters(:, :), parameters_dot(:, :)
    real(dp), allocatable :: vertices(:, :), vertices_dot(:, :)
    real(dp) :: edge_vertices_bar(3, 2), panel_parameters_bar(2, 3)
    integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
    integer, allocatable :: triangles(:, :)
    integer :: basis, panel, status, status_minus, status_plus, vertex
    logical :: all_passed

    all_passed = .true.
    major_radius_dot = 0.013_dp
    minor_radius_dot = -0.009_dp
    xi_dot = 0.071_dp
    eta_dot = -0.053_dp
    point_bar = [0.31_dp, -0.27_dp, 0.19_dp]
    value_bar = [-0.22_dp, 0.17_dp, 0.29_dp]
    divergence_bar = -0.13_dp
    jacobian_bar = 0.37_dp
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 5, 6, vertices, triangles, parameters)
    call build_maxwell_rwg_surface_space( &
        vertices, triangles, edge_vertices, edge_triangles, status)
    call record_condition(status == 0, "torus RWG surface space builds")
    basis = 1
    panel = edge_triangles(1, basis)
    allocate( &
        vertices_dot(size(vertices, 1), size(vertices, 2)), &
        parameters_dot(size(parameters, 1), size(parameters, 2)))
    do vertex = 1, size(vertices, 2)
        vertices_dot(:, vertex) = [ &
            0.007_dp*sin(real(vertex, dp)), &
            -0.005_dp*cos(real(2*vertex, dp)), &
            0.006_dp*sin(real(3*vertex, dp))]
        parameters_dot(:, vertex) = [ &
            0.004_dp*cos(real(vertex + 1, dp)), &
            -0.003_dp*sin(real(2*vertex + 1, dp))]
    end do

    call evaluate_maxwell_torus_curved_rwg_basis_jvp( &
        vertices, triangles, parameters, edge_vertices, edge_triangles, basis, &
        panel, major_radius, minor_radius, xi, eta, vertices_dot, parameters_dot, &
        major_radius_dot, minor_radius_dot, point, value, divergence, jacobian, &
        point_dot, value_dot, divergence_dot, jacobian_dot, status, xi_dot, eta_dot)
    call evaluate_maxwell_torus_curved_rwg_basis( &
        vertices + step*vertices_dot, triangles, parameters + step*parameters_dot, &
        edge_vertices, edge_triangles, basis, panel, &
        major_radius + step*major_radius_dot, minor_radius + step*minor_radius_dot, &
        xi + step*xi_dot, eta + step*eta_dot, point_plus, value_plus, &
        divergence_plus, jacobian_plus, status_plus)
    call evaluate_maxwell_torus_curved_rwg_basis( &
        vertices - step*vertices_dot, triangles, parameters - step*parameters_dot, &
        edge_vertices, edge_triangles, basis, panel, &
        major_radius - step*major_radius_dot, minor_radius - step*minor_radius_dot, &
        xi - step*xi_dot, eta - step*eta_dot, point_minus, value_minus, &
        divergence_minus, jacobian_minus, status_minus)
    call record_condition(status == 0, "torus parent RWG JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed torus parent RWG evaluations succeed")
    call record_condition(maxval(abs(point_dot - &
        (point_plus - point_minus)/(2.0_dp*step))) < 3.0e-8_dp, &
        "torus parent RWG point JVP matches reassembly")
    call record_condition(maxval(abs(value_dot - &
        (value_plus - value_minus)/(2.0_dp*step))) < 3.0e-8_dp, &
        "torus parent RWG value JVP matches reassembly")
    call record_condition(abs(divergence_dot - &
        (divergence_plus - divergence_minus)/(2.0_dp*step)) < 3.0e-8_dp, &
        "torus parent RWG divergence JVP matches reassembly")
    call record_condition(abs(jacobian_dot - &
        (jacobian_plus - jacobian_minus)/(2.0_dp*step)) < 3.0e-8_dp, &
        "torus parent RWG Jacobian JVP matches reassembly")

    call evaluate_maxwell_torus_curved_rwg_basis_vjp( &
        vertices, triangles, parameters, edge_vertices, edge_triangles, basis, &
        panel, major_radius, minor_radius, xi, eta, point_bar, value_bar, &
        divergence_bar, jacobian_bar, edge_vertices_bar, panel_parameters_bar, &
        major_radius_bar, minor_radius_bar, status, xi_bar, eta_bar)
    lhs = dot_product(point_bar, point_dot) + dot_product(value_bar, value_dot) + &
        divergence_bar*divergence_dot + jacobian_bar*jacobian_dot
    rhs = sum(edge_vertices_bar*vertices_dot(:, edge_vertices(:, basis))) + &
        sum(panel_parameters_bar*parameters_dot(:, triangles(:, panel))) + &
        major_radius_bar*major_radius_dot + minor_radius_bar*minor_radius_dot + &
        xi_bar*xi_dot + eta_bar*eta_dot
    call record_condition(status == 0, "torus parent RWG VJP succeeds")
    call record_condition(abs(lhs - rhs) < &
        5.0e-12_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "torus parent RWG products obey the adjoint identity")

    call check_summary("Differentiable exact-curved torus parent RWG basis")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_rwg_basis_ad
