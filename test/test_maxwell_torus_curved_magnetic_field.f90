program test_maxwell_torus_curved_magnetic_field
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_torus_curved_rwg_mass_matrix, &
        evaluate_maxwell_torus_curved_magnetic_field_rwg_3d, &
        evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_jvp, &
        evaluate_maxwell_torus_magnetic_geometry_jvp, &
        evaluate_maxwell_torus_magnetic_geometry_vjp, &
        evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_vjp, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: observation(3) = [0.25_dp, 0.10_dp, 0.15_dp]
    real(dp), allocatable :: vertices(:, :), parameters(:, :), mass(:, :)
    real(dp), allocatable :: vertices_dot(:, :), parameters_dot(:, :)
    integer, allocatable :: triangles(:, :)
    complex(dp), allocatable :: coefficients(:), coefficients_dot(:)
    complex(dp), allocatable :: coefficients_bar(:), doubled(:)
    complex(dp) :: field(3), doubled_field(3), zero_field(3)
    complex(dp) :: field_dot(3), field_bar(3), field_plus(3), field_minus(3)
    complex(dp) :: geometry_dot(3), geometry_plus(3), geometry_minus(3)
    real(dp) :: scaling_error, derivative_error, adjoint_error, geometry_error
    real(dp) :: step
    real(dp) :: major_radius_dot, minor_radius_dot, wave_number_dot
    real(dp) :: observation_dot(3), wave_number
    real(dp), allocatable :: vertices_bar(:, :), parameters_bar(:, :)
    real(dp) :: major_radius_bar, minor_radius_bar, wave_number_bar
    real(dp) :: observation_bar(3), geometry_adjoint_error, lhs, rhs
    complex(dp), allocatable :: coefficients_geometry_bar(:)
    integer :: i, status
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 3, vertices, triangles, parameters)
    call assemble_maxwell_torus_curved_rwg_mass_matrix( &
        vertices, triangles, parameters, major_radius, minor_radius, 2, &
        mass, status)
    call record_condition(status == 0 .and. size(mass, 1) > 0, &
        'curved torus RWG space has a mass matrix')
    allocate( &
        coefficients(size(mass, 1)), coefficients_dot(size(mass, 1)), &
        doubled(size(mass, 1)))
    do i = 1, size(coefficients)
        coefficients(i) = cmplx( &
            cos(0.17_dp*real(i, dp)), sin(0.11_dp*real(i, dp)), dp)
        coefficients_dot(i) = cmplx( &
            0.13_dp*sin(0.07_dp*real(i, dp)), &
            -0.09_dp*cos(0.19_dp*real(i, dp)), dp)
    end do
    doubled = 2.0_dp*coefficients

    call evaluate_maxwell_torus_curved_magnetic_field_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        coefficients, observation, 0.45_dp, 2, field, status)
    call record_condition(status == 0 .and. all(ieee_is_finite(real(field))) &
        .and. all(ieee_is_finite(aimag(field))), &
        'curved torus magnetic field is finite off the surface')
    call evaluate_maxwell_torus_curved_magnetic_field_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, doubled, &
        observation, 0.45_dp, 2, doubled_field, status)
    scaling_error = maxval(abs(doubled_field - 2.0_dp*field))
    call record_condition(status == 0 .and. scaling_error < 5.0e-12_dp, &
        'curved torus magnetic field is linear in RWG coefficients')
    call evaluate_maxwell_torus_curved_magnetic_field_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        cmplx(0.0_dp, 0.0_dp, dp)*coefficients, observation, 0.45_dp, 2, &
        zero_field, status)
    call record_condition(status == 0 .and. maxval(abs(zero_field)) == 0.0_dp, &
        'zero curved torus current has zero magnetic field')

    step = 1.0e-6_dp
    call evaluate_maxwell_torus_curved_magnetic_field_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        coefficients + step*coefficients_dot, observation, 0.45_dp, 2, &
        field_plus, status)
    call evaluate_maxwell_torus_curved_magnetic_field_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        coefficients - step*coefficients_dot, observation, 0.45_dp, 2, &
        field_minus, status)
    call evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_jvp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        coefficients_dot, observation, 0.45_dp, 2, field_dot, status)
    derivative_error = maxval(abs(field_dot - &
        (field_plus - field_minus)/(2.0_dp*step)))
    call record_condition(status == 0 .and. derivative_error < 2.0e-8_dp, &
        'curved torus magnetic-field JVP matches finite differences')

    field_bar = [ &
        cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.1_dp, 0.6_dp, dp)]
    call evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_vjp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        observation, 0.45_dp, 2, field_bar, coefficients_bar, status)
    adjoint_error = abs(real(sum(conjg(field_bar)*field_dot)) - &
        real(sum(conjg(coefficients_bar)*coefficients_dot)))
    call record_condition(status == 0 .and. adjoint_error < 2.0e-11_dp, &
        'curved torus magnetic-field VJP satisfies the dot-product identity')

    allocate(vertices_dot(size(vertices, 1), size(vertices, 2)))
    allocate(parameters_dot(size(parameters, 1), size(parameters, 2)))
    vertices_dot = 0.0_dp
    parameters_dot = 0.0_dp
    do i = 1, size(vertices, 2)
        vertices_dot(:, i) = [0.003_dp*i, -0.002_dp*i, 0.001_dp*i]
        parameters_dot(:, i) = [0.0007_dp*i, -0.0004_dp*i]
    end do
    major_radius_dot = 0.002_dp
    minor_radius_dot = -0.001_dp
    observation_dot = [0.004_dp, -0.003_dp, 0.002_dp]
    wave_number = 0.45_dp
    wave_number_dot = 0.017_dp
    call evaluate_maxwell_torus_magnetic_geometry_jvp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        coefficients, observation, wave_number, 2, vertices_dot, &
        parameters_dot, major_radius_dot, minor_radius_dot, coefficients_dot, &
        observation_dot, wave_number_dot, geometry_dot, status)
    call evaluate_maxwell_torus_curved_magnetic_field_rwg_3d( &
        vertices + step*vertices_dot, triangles, parameters + step*parameters_dot, &
        major_radius + step*major_radius_dot, minor_radius + step*minor_radius_dot, &
        coefficients + step*coefficients_dot, observation + step*observation_dot, &
        wave_number + step*wave_number_dot, 2, geometry_plus, status)
    call evaluate_maxwell_torus_curved_magnetic_field_rwg_3d( &
        vertices - step*vertices_dot, triangles, parameters - step*parameters_dot, &
        major_radius - step*major_radius_dot, minor_radius - step*minor_radius_dot, &
        coefficients - step*coefficients_dot, observation - step*observation_dot, &
        wave_number - step*wave_number_dot, 2, geometry_minus, status)
    geometry_error = maxval(abs(geometry_dot - &
        (geometry_plus - geometry_minus)/(2.0_dp*step)))
    call record_condition(status == 0 .and. geometry_error < 2.0e-7_dp, &
        'curved torus magnetic-field geometry JVP matches reassembly')
    allocate( &
        vertices_bar(size(vertices, 1), size(vertices, 2)), &
        parameters_bar(size(parameters, 1), size(parameters, 2)))
    call evaluate_maxwell_torus_magnetic_geometry_vjp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        coefficients, observation, wave_number, 2, field_bar, vertices_bar, &
        parameters_bar, major_radius_bar, minor_radius_bar, &
        coefficients_geometry_bar, observation_bar, wave_number_bar, status)
    lhs = real(sum(conjg(field_bar)*geometry_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + sum(parameters_bar*parameters_dot) + &
        major_radius_bar*major_radius_dot + minor_radius_bar*minor_radius_dot + &
        real(sum(conjg(coefficients_geometry_bar)*coefficients_dot), dp) + &
        dot_product(observation_bar, observation_dot) + wave_number_bar* &
        wave_number_dot
    geometry_adjoint_error = abs(lhs - rhs)
    call record_condition(status == 0, &
        'curved torus magnetic-field geometry VJP succeeds')
    call record_condition(geometry_adjoint_error < &
        5.0e-9_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        'curved torus magnetic-field geometry VJP satisfies the adjoint identity')

    call check_summary('Exact-curved torus RWG magnetic field')
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_magnetic_field
