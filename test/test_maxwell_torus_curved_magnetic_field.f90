program test_maxwell_torus_curved_magnetic_field
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_torus_curved_rwg_mass_matrix, &
        evaluate_maxwell_torus_curved_magnetic_field_rwg_3d, &
        evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_jvp, &
        evaluate_maxwell_torus_curved_magnetic_field_rwg_3d_vjp, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: observation(3) = [0.25_dp, 0.10_dp, 0.15_dp]
    real(dp), allocatable :: vertices(:, :), parameters(:, :), mass(:, :)
    integer, allocatable :: triangles(:, :)
    complex(dp), allocatable :: coefficients(:), coefficients_dot(:)
    complex(dp), allocatable :: coefficients_bar(:), doubled(:)
    complex(dp) :: field(3), doubled_field(3), zero_field(3)
    complex(dp) :: field_dot(3), field_bar(3), field_plus(3), field_minus(3)
    real(dp) :: scaling_error, derivative_error, adjoint_error, step
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
