program test_torus_curved_helmholtz_dtn_geometry_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_helmholtz_torus_curved_dtn_3d, &
        assemble_helmholtz_torus_curved_dtn_3d_geometry_jvp, &
        assemble_helmholtz_torus_curved_dtn_3d_geometry_vjp, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.35_dp
    real(dp), parameter :: step = 2.0e-6_dp
    integer, parameter :: quadrature_degree = 3
    integer, parameter :: count_theta = 4, count_phi = 6
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), parameters_dot(:, :)
    real(dp), allocatable :: parameters_bar(:, :), vertices(:, :)
    complex(dp), allocatable :: single_layer(:, :), double_layer(:, :)
    complex(dp), allocatable :: single_layer_dot(:, :), double_layer_dot(:, :)
    complex(dp), allocatable :: single_layer_plus(:, :), double_layer_plus(:, :)
    complex(dp), allocatable :: single_layer_minus(:, :), double_layer_minus(:, :)
    complex(dp), allocatable :: single_layer_bar(:, :), double_layer_bar(:, :)
    real(dp), allocatable :: mass(:, :), mass_dot(:, :), mass_plus(:, :)
    real(dp), allocatable :: mass_minus(:, :), mass_bar(:, :)
    real(dp) :: major_radius_dot, minor_radius_dot, wave_number_dot
    real(dp) :: major_radius_bar, minor_radius_bar, wave_number_bar
    real(dp) :: relative_difference, lhs, rhs
    integer :: status, status_minus, status_plus, vertex
    logical :: all_passed

    all_passed = .true.
    major_radius_dot = 0.015_dp
    minor_radius_dot = -0.01_dp
    wave_number_dot = 0.02_dp
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, count_theta, count_phi, vertices, &
        triangles, parameters)
    allocate(parameters_dot(size(parameters, 1), size(parameters, 2)))
    allocate(parameters_bar(size(parameters, 1), size(parameters, 2)))
    do vertex = 1, size(parameters, 2)
        parameters_dot(:, vertex) = [ &
            0.006_dp*cos(parameters(1, vertex)), &
            -0.004_dp*sin(parameters(2, vertex))]
    end do

    call assemble_helmholtz_torus_curved_dtn_3d_geometry_jvp( &
        parameters, triangles, major_radius, minor_radius, wave_number, &
        quadrature_degree, parameters_dot, major_radius_dot, &
        minor_radius_dot, wave_number_dot, single_layer_dot, double_layer_dot, &
        mass_dot, status)
    call assemble_helmholtz_torus_curved_dtn_3d( &
        parameters + step*parameters_dot, triangles, &
        major_radius + step*major_radius_dot, minor_radius + step*minor_radius_dot, &
        wave_number + step*wave_number_dot, quadrature_degree, &
        single_layer_plus, double_layer_plus, mass_plus, status_plus)
    call assemble_helmholtz_torus_curved_dtn_3d( &
        parameters - step*parameters_dot, triangles, &
        major_radius - step*major_radius_dot, minor_radius - step*minor_radius_dot, &
        wave_number - step*wave_number_dot, quadrature_degree, &
        single_layer_minus, double_layer_minus, mass_minus, status_minus)
    relative_difference = max( &
        maxval(abs(single_layer_dot - (single_layer_plus - single_layer_minus)/ &
        (2.0_dp*step))), &
        max(maxval(abs(double_layer_dot - &
        (double_layer_plus - double_layer_minus)/(2.0_dp*step))), &
        maxval(abs(mass_dot - (mass_plus - mass_minus)/(2.0_dp*step)))))
    call record_condition(status == 0 .and. status_plus == 0 .and. &
        status_minus == 0, "toroidal Helmholtz DtN geometry JVP succeeds")
    call record_condition(relative_difference < 3.0e-7_dp, &
        "toroidal Helmholtz DtN geometry JVP matches finite differences")

    allocate(single_layer(size(single_layer_dot, 1), size(single_layer_dot, 2)), &
        double_layer(size(double_layer_dot, 1), size(double_layer_dot, 2)), &
        mass(size(mass_dot, 1), size(mass_dot, 2)))
    call assemble_helmholtz_torus_curved_dtn_3d( &
        parameters, triangles, major_radius, minor_radius, wave_number, &
        quadrature_degree, single_layer, double_layer, mass, status)
    allocate(single_layer_bar(size(single_layer, 1), size(single_layer, 2)), &
        double_layer_bar(size(double_layer, 1), size(double_layer, 2)), &
        mass_bar(size(mass, 1), size(mass, 2)))
    single_layer_bar = cmplx(0.1_dp, -0.07_dp, dp)
    double_layer_bar = cmplx(0.2_dp, 0.04_dp, dp)
    mass_bar = 0.3_dp
    call assemble_helmholtz_torus_curved_dtn_3d_geometry_vjp( &
        parameters, triangles, major_radius, minor_radius, wave_number, &
        quadrature_degree, single_layer_bar, double_layer_bar, mass_bar, &
        parameters_bar, major_radius_bar, minor_radius_bar, wave_number_bar, &
        status)
    lhs = sum(real(conjg(single_layer_bar)*single_layer_dot, dp)) + &
        sum(real(conjg(double_layer_bar)*double_layer_dot, dp)) + &
        sum(mass_bar*mass_dot)
    rhs = sum(parameters_bar*parameters_dot) + &
        major_radius_bar*major_radius_dot + minor_radius_bar*minor_radius_dot + &
        wave_number_bar*wave_number_dot
    call record_condition(status == 0, &
        "toroidal Helmholtz DtN geometry VJP succeeds")
    call record_condition(abs(lhs - rhs) < 3.0e-9_dp*max(1.0_dp, &
        abs(lhs), abs(rhs)), &
        "toroidal Helmholtz DtN geometry products obey the adjoint identity")

    call check_summary("Torus curved Helmholtz DtN geometry derivatives")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_curved_helmholtz_dtn_geometry_ad
