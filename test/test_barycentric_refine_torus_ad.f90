program test_barycentric_refine_torus_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        barycentric_refine_torus_surface_mesh, &
        barycentric_refine_torus_surface_mesh_jvp, &
        barycentric_refine_torus_surface_mesh_vjp, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.3_dp, minor_radius = 0.57_dp
    real(dp), parameter :: step = 2.0e-6_dp
    real(dp) :: major_radius_bar, minor_radius_bar
    real(dp) :: major_radius_dot, minor_radius_dot
    real(dp) :: lhs, rhs
    real(dp), allocatable :: parameters(:, :), parameters_bar(:, :)
    real(dp), allocatable :: parameters_dot(:, :), vertices(:, :), vertices_bar(:, :)
    real(dp), allocatable :: vertices_dot(:, :)
    real(dp), allocatable :: refined_parameters(:, :), refined_parameters_bar(:, :)
    real(dp), allocatable :: refined_parameters_dot(:, :)
    real(dp), allocatable :: refined_vertices(:, :), refined_vertices_bar(:, :)
    real(dp), allocatable :: refined_vertices_dot(:, :)
    real(dp), allocatable :: refined_vertices_minus(:, :), refined_vertices_plus(:, :)
    real(dp), allocatable :: refined_parameters_minus(:, :), refined_parameters_plus(:, :)
    integer, allocatable :: triangles(:, :), refined_triangles(:, :)
    integer, allocatable :: refined_triangles_jvp(:, :)
    integer :: status, status_minus, status_plus, vertex
    logical :: all_passed

    all_passed = .true.
    major_radius_dot = 0.023_dp
    minor_radius_dot = -0.011_dp
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 4, 5, vertices, triangles, parameters)
    allocate( &
        vertices_dot(size(vertices, 1), size(vertices, 2)), &
        parameters_dot(size(parameters, 1), size(parameters, 2)), &
        vertices_bar(size(vertices, 1), size(vertices, 2)), &
        parameters_bar(size(parameters, 1), size(parameters, 2)))
    do vertex = 1, size(vertices, 2)
        vertices_dot(:, vertex) = [ &
            0.006_dp*sin(real(vertex, dp)), &
            -0.004_dp*cos(real(2*vertex, dp)), &
            0.005_dp*sin(real(3*vertex, dp))]
        parameters_dot(:, vertex) = [ &
            0.005_dp*cos(real(vertex + 1, dp)), &
            -0.003_dp*sin(real(2*vertex + 1, dp))]
    end do
    call barycentric_refine_torus_surface_mesh( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        refined_vertices, refined_triangles, refined_parameters, status)
    call barycentric_refine_torus_surface_mesh_jvp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        vertices_dot, parameters_dot, major_radius_dot, minor_radius_dot, &
        refined_vertices_dot, refined_triangles_jvp, refined_parameters_dot, &
        status)
    call barycentric_refine_torus_surface_mesh( &
        vertices + step*vertices_dot, triangles, parameters + step*parameters_dot, &
        major_radius + step*major_radius_dot, minor_radius + step*minor_radius_dot, &
        refined_vertices_plus, refined_triangles, refined_parameters_plus, &
        status_plus)
    call barycentric_refine_torus_surface_mesh( &
        vertices - step*vertices_dot, triangles, parameters - step*parameters_dot, &
        major_radius - step*major_radius_dot, minor_radius - step*minor_radius_dot, &
        refined_vertices_minus, refined_triangles, refined_parameters_minus, &
        status_minus)
    call record_condition(status == 0, "torus barycentric refinement succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed torus barycentric refinements succeed")
    call record_condition(all(refined_triangles_jvp == refined_triangles), &
        "refinement JVP preserves the discrete topology")
    call record_condition(maxval(abs(refined_vertices_dot - &
        (refined_vertices_plus - refined_vertices_minus)/(2.0_dp*step))) < 4.0e-8_dp, &
        "torus refinement vertex JVP matches reassembly")
    call record_condition(maxval(abs(refined_parameters_dot - &
        (refined_parameters_plus - refined_parameters_minus)/(2.0_dp*step))) < &
        4.0e-8_dp, "torus refinement parameter JVP matches reassembly")

    allocate( &
        refined_vertices_bar(size(refined_vertices, 1), size(refined_vertices, 2)), &
        refined_parameters_bar(size(refined_parameters, 1), size(refined_parameters, 2)))
    do vertex = 1, size(refined_vertices_bar, 2)
        refined_vertices_bar(:, vertex) = [ &
            0.007_dp*sin(real(vertex + 1, dp)), &
            -0.005_dp*cos(real(2*vertex + 1, dp)), &
            0.006_dp*sin(real(3*vertex + 1, dp))]
        refined_parameters_bar(:, vertex) = [ &
            0.004_dp*cos(real(vertex, dp)), &
            -0.003_dp*sin(real(2*vertex, dp))]
    end do
    call barycentric_refine_torus_surface_mesh_vjp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        refined_vertices_bar, refined_parameters_bar, vertices_bar, &
        parameters_bar, major_radius_bar, minor_radius_bar, status)
    lhs = sum(refined_vertices_bar*refined_vertices_dot) + &
        sum(refined_parameters_bar*refined_parameters_dot)
    rhs = sum(vertices_bar*vertices_dot) + sum(parameters_bar*parameters_dot) + &
        major_radius_bar*major_radius_dot + minor_radius_bar*minor_radius_dot
    call record_condition(status == 0, "torus refinement VJP succeeds")
    call record_condition(abs(lhs - rhs) < 5.0e-12_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "torus refinement products obey the adjoint identity")

    call check_summary("Differentiable torus barycentric refinement")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_barycentric_refine_torus_ad
