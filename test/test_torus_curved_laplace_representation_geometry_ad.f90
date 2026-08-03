program test_torus_curved_laplace_representation_geometry_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        evaluate_laplace_representation_torus_curved_3d, &
        evaluate_laplace_representation_torus_curved_3d_geometry_jvp, &
        evaluate_laplace_representation_torus_curved_3d_geometry_vjp
    use fortfem_core, only: &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: step = 2.0e-6_dp
    integer, parameter :: quadrature_degree = 4
    integer, parameter :: count_theta = 6, count_phi = 8
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), parameters_dot(:, :)
    real(dp), allocatable :: parameters_bar(:, :), vertices(:, :)
    real(dp), allocatable :: trace(:), flux(:)
    real(dp) :: target(3), target_dot(3), target_bar(3)
    real(dp) :: value_dot, value_bar, value_plus, value_minus
    real(dp) :: finite_difference, lhs, rhs
    real(dp) :: major_radius_dot, minor_radius_dot
    real(dp) :: major_radius_bar, minor_radius_bar
    integer :: status, status_minus, status_plus, vertex
    logical :: all_passed

    all_passed = .true.
    target = [0.25_dp, 0.17_dp, 0.13_dp]
    target_dot = [-0.03_dp, 0.04_dp, 0.02_dp]
    major_radius_dot = 0.015_dp
    minor_radius_dot = -0.01_dp
    value_bar = -1.7_dp
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, count_theta, count_phi, vertices, &
        triangles, parameters)
    allocate(parameters_dot(size(parameters, 1), size(parameters, 2)))
    allocate(parameters_bar(size(parameters, 1), size(parameters, 2)))
    allocate(trace(size(vertices, 2)), flux(size(vertices, 2)))
    do vertex = 1, size(vertices, 2)
        parameters_dot(:, vertex) = [ &
            0.006_dp*cos(parameters(1, vertex)), &
            -0.004_dp*sin(parameters(2, vertex))]
        trace(vertex) = 0.7_dp + 0.11_dp*parameters(1, vertex) - &
            0.08_dp*parameters(2, vertex)
        flux(vertex) = -0.2_dp + 0.06_dp*cos(parameters(1, vertex)) + &
            0.05_dp*sin(parameters(2, vertex))
    end do

    call evaluate_laplace_representation_torus_curved_3d_geometry_jvp( &
        parameters, triangles, major_radius, minor_radius, trace, flux, &
        target, quadrature_degree, parameters_dot, major_radius_dot, &
        minor_radius_dot, target_dot, value_dot, status)
    call evaluate_laplace_representation_torus_curved_3d( &
        parameters + step*parameters_dot, triangles, &
        major_radius + step*major_radius_dot, minor_radius + step*minor_radius_dot, &
        trace, flux, target + step*target_dot, quadrature_degree, value_plus, &
        status_plus)
    call evaluate_laplace_representation_torus_curved_3d( &
        parameters - step*parameters_dot, triangles, &
        major_radius - step*major_radius_dot, minor_radius - step*minor_radius_dot, &
        trace, flux, target - step*target_dot, quadrature_degree, value_minus, &
        status_minus)
    finite_difference = (value_plus - value_minus)/(2.0_dp*step)
    call record_condition(status == 0 .and. status_plus == 0 .and. &
        status_minus == 0, "toroidal Laplace geometry JVP succeeds")
    call record_condition(abs(value_dot - finite_difference) < 2.0e-8_dp, &
        "toroidal Laplace geometry JVP matches finite differences")

    call evaluate_laplace_representation_torus_curved_3d_geometry_vjp( &
        parameters, triangles, major_radius, minor_radius, trace, flux, &
        target, quadrature_degree, value_bar, parameters_bar, &
        major_radius_bar, minor_radius_bar, target_bar, status)
    lhs = value_bar*value_dot
    rhs = sum(parameters_bar*parameters_dot) + &
        major_radius_bar*major_radius_dot + minor_radius_bar*minor_radius_dot + &
        dot_product(target_bar, target_dot)
    call record_condition(status == 0, "toroidal Laplace geometry VJP succeeds")
    call record_condition(abs(lhs - rhs) < 2.0e-10_dp*max(1.0_dp, &
        abs(lhs), abs(rhs)), &
        "toroidal Laplace geometry products obey the adjoint identity")

    call check_summary("Torus curved Laplace representation geometry derivatives")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_curved_laplace_representation_geometry_ad
