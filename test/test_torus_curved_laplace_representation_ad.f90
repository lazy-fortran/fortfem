program test_torus_curved_laplace_representation_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_laplace_representation_torus_curved_3d, &
        evaluate_laplace_representation_torus_curved_3d_jvp, &
        evaluate_laplace_representation_torus_curved_3d_vjp, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: step = 2.0e-6_dp
    integer, parameter :: quadrature_degree = 4
    integer, parameter :: count_theta = 6, count_phi = 8
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    real(dp), allocatable :: trace(:), trace_dot(:), trace_bar(:)
    real(dp), allocatable :: flux(:), flux_dot(:), flux_bar(:)
    real(dp), allocatable :: panel_flux(:), panel_flux_dot(:)
    real(dp), allocatable :: panel_flux_bar(:)
    real(dp) :: target(3), target_dot(3), target_bar(3)
    real(dp) :: value, value_dot, value_bar, value_plus, value_minus
    real(dp) :: panel_value_dot, panel_value_plus, panel_value_minus
    real(dp) :: finite_difference, panel_finite_difference
    real(dp) :: lhs, rhs, panel_lhs, panel_rhs
    integer :: status, status_minus, status_plus, vertex, element
    logical :: all_passed

    all_passed = .true.
    target = [0.25_dp, 0.17_dp, 0.13_dp]
    target_dot = [-0.03_dp, 0.04_dp, 0.02_dp]
    value_bar = -1.7_dp
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, count_theta, count_phi, vertices, &
        triangles, parameters)
    allocate(trace(size(vertices, 2)), trace_dot(size(vertices, 2)))
    allocate(trace_bar(size(vertices, 2)))
    allocate(flux(size(vertices, 2)), flux_dot(size(vertices, 2)))
    allocate(flux_bar(size(vertices, 2)))
    do vertex = 1, size(vertices, 2)
        trace(vertex) = 0.7_dp + 0.11_dp*parameters(1, vertex) - &
            0.08_dp*parameters(2, vertex)
        trace_dot(vertex) = 0.03_dp*cos(parameters(1, vertex)) + &
            0.02_dp*sin(parameters(2, vertex))
        flux(vertex) = -0.2_dp + 0.06_dp*cos(parameters(1, vertex)) + &
            0.05_dp*sin(parameters(2, vertex))
        flux_dot(vertex) = -0.04_dp*sin(parameters(1, vertex)) + &
            0.01_dp*cos(parameters(2, vertex))
    end do

    call evaluate_laplace_representation_torus_curved_3d_jvp( &
        parameters, triangles, major_radius, minor_radius, trace, flux, &
        target, quadrature_degree, trace_dot, flux_dot, target_dot, &
        value_dot, status)
    call evaluate_laplace_representation_torus_curved_3d( &
        parameters, triangles, major_radius, minor_radius, trace + &
        step*trace_dot, flux + step*flux_dot, target + step*target_dot, &
        quadrature_degree, value_plus, status_plus)
    call evaluate_laplace_representation_torus_curved_3d( &
        parameters, triangles, major_radius, minor_radius, trace - &
        step*trace_dot, flux - step*flux_dot, target - step*target_dot, &
        quadrature_degree, value_minus, status_minus)
    finite_difference = (value_plus - value_minus)/(2.0_dp*step)
    call record_condition(status == 0 .and. status_plus == 0 .and. &
        status_minus == 0, "toroidal Laplace representation JVP succeeds")
    call record_condition(abs(value_dot - finite_difference) < 2.0e-8_dp, &
        "toroidal Laplace representation JVP matches finite differences")

    call evaluate_laplace_representation_torus_curved_3d_vjp( &
        parameters, triangles, major_radius, minor_radius, trace, flux, &
        target, quadrature_degree, value_bar, trace_bar, flux_bar, &
        target_bar, status)
    lhs = value_bar*value_dot
    rhs = dot_product(trace_bar, trace_dot) + dot_product(flux_bar, flux_dot) &
        + dot_product(target_bar, target_dot)
    call record_condition(status == 0, &
        "toroidal Laplace representation VJP succeeds")
    call record_condition(abs(lhs - rhs) < 2.0e-10_dp*max(1.0_dp, &
        abs(lhs), abs(rhs)), &
        "toroidal Laplace representation products obey the adjoint identity")

    allocate(panel_flux(size(triangles, 2)), &
        panel_flux_dot(size(triangles, 2)), panel_flux_bar(size(triangles, 2)))
    do element = 1, size(triangles, 2)
        panel_flux(element) = 0.1_dp + 0.003_dp*real(element, dp)
        panel_flux_dot(element) = -0.02_dp + 0.001_dp*real(element, dp)
    end do
    call evaluate_laplace_representation_torus_curved_3d_jvp( &
        parameters, triangles, major_radius, minor_radius, trace, panel_flux, &
        target, quadrature_degree, trace_dot, panel_flux_dot, target_dot, &
        panel_value_dot, status)
    call evaluate_laplace_representation_torus_curved_3d( &
        parameters, triangles, major_radius, minor_radius, trace + &
        step*trace_dot, panel_flux + step*panel_flux_dot, &
        target + step*target_dot, quadrature_degree, panel_value_plus, &
        status_plus)
    call evaluate_laplace_representation_torus_curved_3d( &
        parameters, triangles, major_radius, minor_radius, trace - &
        step*trace_dot, panel_flux - step*panel_flux_dot, &
        target - step*target_dot, quadrature_degree, panel_value_minus, &
        status_minus)
    panel_finite_difference = (panel_value_plus - panel_value_minus)/ &
        (2.0_dp*step)
    call record_condition(status == 0 .and. status_plus == 0 .and. &
        status_minus == 0, "panel-flux Laplace representation JVP succeeds")
    call record_condition(abs(panel_value_dot - panel_finite_difference) < &
        2.0e-8_dp, "panel-flux Laplace representation JVP matches differences")
    call evaluate_laplace_representation_torus_curved_3d_vjp( &
        parameters, triangles, major_radius, minor_radius, trace, panel_flux, &
        target, quadrature_degree, value_bar, trace_bar, panel_flux_bar, &
        target_bar, status)
    panel_lhs = value_bar*panel_value_dot
    panel_rhs = dot_product(trace_bar, trace_dot) + &
        dot_product(panel_flux_bar, panel_flux_dot) + &
        dot_product(target_bar, target_dot)
    call record_condition(status == 0, &
        "panel-flux Laplace representation VJP succeeds")
    call record_condition(abs(panel_lhs - panel_rhs) < &
        2.0e-10_dp*max(1.0_dp, abs(panel_lhs), abs(panel_rhs)), &
        "panel-flux representation products obey the adjoint identity")

    call check_summary("Torus curved Laplace representation derivatives")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_curved_laplace_representation_ad
