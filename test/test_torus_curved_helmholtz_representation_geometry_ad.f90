program test_torus_curved_helmholtz_representation_geometry_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_helmholtz_representation_torus_curved_3d, &
        evaluate_helmholtz_representation_torus_curved_3d_geometry_jvp, &
        evaluate_helmholtz_representation_torus_curved_3d_geometry_vjp, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.7_dp, step = 2.0e-6_dp
    integer, parameter :: quadrature_degree = 4
    integer, parameter :: count_theta = 6, count_phi = 8
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), parameters_dot(:, :)
    real(dp), allocatable :: parameters_bar(:, :), vertices(:, :)
    complex(dp), allocatable :: trace(:), flux(:)
    real(dp) :: target(3), target_dot(3), target_bar(3)
    real(dp) :: major_radius_dot, minor_radius_dot
    real(dp) :: major_radius_bar, minor_radius_bar
    complex(dp) :: value_dot, value_bar, value_plus, value_minus
    real(dp) :: finite_difference, lhs, rhs
    integer :: status, status_minus, status_plus, vertex
    logical :: all_passed

    all_passed = .true.
    target = [0.25_dp, 0.17_dp, 0.13_dp]
    target_dot = [-0.03_dp, 0.04_dp, 0.02_dp]
    major_radius_dot = 0.015_dp
    minor_radius_dot = -0.01_dp
    value_bar = cmplx(-0.6_dp, 0.35_dp, dp)
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
        trace(vertex) = cmplx( &
            0.7_dp + 0.11_dp*parameters(1, vertex), &
            -0.08_dp*parameters(2, vertex), dp)
        flux(vertex) = cmplx( &
            -0.2_dp + 0.06_dp*cos(parameters(1, vertex)), &
            0.05_dp*sin(parameters(2, vertex)), dp)
    end do

    call evaluate_helmholtz_representation_torus_curved_3d_geometry_jvp( &
        parameters, triangles, major_radius, minor_radius, trace, flux, &
        target, wave_number, quadrature_degree, parameters_dot, &
        major_radius_dot, minor_radius_dot, target_dot, value_dot, status)
    call evaluate_helmholtz_representation_torus_curved_3d( &
        parameters + step*parameters_dot, triangles, &
        major_radius + step*major_radius_dot, minor_radius + step*minor_radius_dot, &
        trace, flux, target + step*target_dot, wave_number, quadrature_degree, &
        value_plus, status_plus)
    call evaluate_helmholtz_representation_torus_curved_3d( &
        parameters - step*parameters_dot, triangles, &
        major_radius - step*major_radius_dot, minor_radius - step*minor_radius_dot, &
        trace, flux, target - step*target_dot, wave_number, quadrature_degree, &
        value_minus, status_minus)
    finite_difference = abs( &
        (value_plus - value_minus)/(2.0_dp*step) - value_dot)
    call record_condition(status == 0 .and. status_plus == 0 .and. &
        status_minus == 0, "toroidal Helmholtz geometry JVP succeeds")
    call record_condition(finite_difference < 3.0e-8_dp, &
        "toroidal Helmholtz geometry JVP matches finite differences")

    call evaluate_helmholtz_representation_torus_curved_3d_geometry_vjp( &
        parameters, triangles, major_radius, minor_radius, trace, flux, &
        target, wave_number, quadrature_degree, value_bar, parameters_bar, &
        major_radius_bar, minor_radius_bar, target_bar, status)
    lhs = real(conjg(value_bar)*value_dot, dp)
    rhs = sum(parameters_bar*parameters_dot) + &
        major_radius_bar*major_radius_dot + minor_radius_bar*minor_radius_dot + &
        dot_product(target_bar, target_dot)
    call record_condition(status == 0, "toroidal Helmholtz geometry VJP succeeds")
    call record_condition(abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, &
        abs(lhs), abs(rhs)), &
        "toroidal Helmholtz geometry products obey the real adjoint identity")

    call check_summary("Torus curved Helmholtz representation geometry derivatives")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_curved_helmholtz_representation_geometry_ad
