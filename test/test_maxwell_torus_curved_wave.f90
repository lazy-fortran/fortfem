program test_maxwell_torus_curved_wave
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d, &
        assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d_jvp, &
        assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d_vjp, &
        evaluate_maxwell_torus_curved_far_field_rwg_3d, &
        evaluate_maxwell_torus_curved_far_field_rwg_3d_jvp, &
        evaluate_maxwell_torus_curved_far_field_rwg_3d_vjp, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: impedance = 1.7_dp, wave_number = 0.45_dp
    complex(dp), allocatable :: coefficients(:), right_hand_side(:)
    complex(dp), allocatable :: scaled_right_hand_side(:)
    complex(dp), allocatable :: right_hand_side_dot(:), right_hand_side_bar(:)
    complex(dp), allocatable :: right_hand_side_plus(:), right_hand_side_minus(:)
    complex(dp), allocatable :: polarization_bar(:)
    complex(dp), allocatable :: coefficients_dot(:), coefficients_bar(:)
    complex(dp) :: far_field_dot(3), far_field_plus(3), far_field_minus(3)
    complex(dp) :: far_field_bar(3)
    complex(dp) :: far_field(3), polarization(3), received, transmitted
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    real(dp), allocatable :: scaled_vertices(:, :)
    real(dp), allocatable :: vertices_dot(:, :), parameters_dot(:, :)
    real(dp), allocatable :: vertices_bar(:, :), parameters_bar(:, :)
    real(dp) :: direction(3), relative_scaling_error
    real(dp) :: direction_dot(3), direction_axis(3), major_radius_dot
    real(dp) :: minor_radius_dot, major_radius_bar, minor_radius_bar
    real(dp) :: direction_bar(3), wave_number_bar, wave_number_dot
    real(dp) :: impedance_dot, impedance_bar
    real(dp) :: step, jvp_error, lhs, rhs, adjoint_error
    complex(dp) :: polarization_dot(3)
    integer, parameter :: derivative_quadrature_degree = 4
    integer :: basis, status
    logical :: all_passed

    all_passed = .true.
    direction = [1.0_dp, 2.0_dp, -1.0_dp]/sqrt(6.0_dp)
    polarization = cmplx([2.0_dp, -1.0_dp, 0.0_dp]/sqrt(5.0_dp), 0.0_dp, dp)
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 5, 6, vertices, triangles, parameters)
    call assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        -direction, polarization, wave_number, 12, right_hand_side, status)
    call record_condition(status == 0, &
        "exact-torus plane-wave RWG trace load assembles")
    allocate(coefficients(size(right_hand_side)))
    do basis = 1, size(coefficients)
        coefficients(basis) = cmplx( &
            cos(real(2*basis, dp)), sin(real(3*basis, dp)), dp)
    end do
    call evaluate_maxwell_torus_curved_far_field_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        coefficients, direction, wave_number, impedance, 12, far_field, status)
    received = sum(polarization*far_field)
    transmitted = -cmplx( &
        0.0_dp, wave_number*impedance/(4.0_dp*acos(-1.0_dp)), dp)* &
        sum(coefficients*right_hand_side)
    call record_condition(status == 0 .and. abs(received - transmitted) < &
        2.0e-13_dp*max(1.0_dp, abs(transmitted)), &
        "exact-torus RWG radiation obeys receive-transmit reciprocity")

    call assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        -direction, polarization, 0.0_dp, 12, right_hand_side, status)
    scaled_vertices = 3.0_dp*vertices
    call assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d( &
        scaled_vertices, triangles, parameters, 3.0_dp*major_radius, &
        3.0_dp*minor_radius, -direction, polarization, 0.0_dp, 12, &
        scaled_right_hand_side, status)
    relative_scaling_error = &
        maxval(abs(scaled_right_hand_side - 9.0_dp*right_hand_side))/ &
        maxval(abs(scaled_right_hand_side))
    call record_condition(status == 0 .and. relative_scaling_error < &
        5.0e-14_dp, "constant-field trace load obeys analytical area scaling")

    allocate( &
        vertices_dot(size(vertices, 1), size(vertices, 2)), &
        parameters_dot(size(parameters, 1), size(parameters, 2)))
    do basis = 1, size(vertices_dot, 2)
        vertices_dot(:, basis) = [ &
            0.011_dp*sin(real(2*basis, dp)), &
            -0.008_dp*cos(real(3*basis, dp)), &
            0.009_dp*sin(real(5*basis, dp))]
        parameters_dot(:, basis) = [ &
            0.006_dp*cos(real(basis + 1, dp)), &
            -0.004_dp*sin(real(2*basis + 1, dp))]
    end do
    direction_axis = [0.3_dp, -0.4_dp, 0.5_dp]
    direction_axis = direction_axis - &
        direction*dot_product(direction, direction_axis)
    direction_axis = direction_axis/norm2(direction_axis)
    direction_dot = real_cross_product(direction_axis, direction)
    polarization_dot = cmplx(real_cross_product( &
        direction_axis, real(polarization, dp)), 0.0_dp, dp)
    major_radius_dot = 0.017_dp
    minor_radius_dot = -0.009_dp
    wave_number_dot = 0.031_dp
    step = 2.0e-6_dp
    call assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d_jvp( &
        vertices, triangles, parameters, major_radius, minor_radius, direction, &
        polarization, wave_number, derivative_quadrature_degree, &
        vertices_dot, parameters_dot, major_radius_dot, minor_radius_dot, &
        direction_dot, polarization_dot, wave_number_dot, right_hand_side, &
        right_hand_side_dot, status)
    call assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d( &
        vertices + step*vertices_dot, triangles, parameters + step*parameters_dot, &
        major_radius + step*major_radius_dot, minor_radius + step*minor_radius_dot, &
        rotate_vector(direction, direction_axis, step), polarization + &
        step*polarization_dot, wave_number + step*wave_number_dot, &
        derivative_quadrature_degree, right_hand_side_plus, status)
    call assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d( &
        vertices - step*vertices_dot, triangles, parameters - step*parameters_dot, &
        major_radius - step*major_radius_dot, minor_radius - step*minor_radius_dot, &
        rotate_vector(direction, direction_axis, -step), polarization - &
        step*polarization_dot, wave_number - step*wave_number_dot, &
        derivative_quadrature_degree, right_hand_side_minus, status)
    jvp_error = maxval(abs(right_hand_side_dot - &
        (right_hand_side_plus - right_hand_side_minus)/(2.0_dp*step)))
    call record_condition(status == 0 .and. jvp_error < 2.0e-7_dp, &
        "exact-torus plane-wave RHS geometry/data JVP matches reassembly")

    allocate( &
        right_hand_side_bar(size(right_hand_side)), &
        vertices_bar(size(vertices, 1), size(vertices, 2)), &
        parameters_bar(size(parameters, 1), size(parameters, 2)))
    do basis = 1, size(right_hand_side_bar)
        right_hand_side_bar(basis) = cmplx( &
            sin(real(2*basis, dp)), cos(real(3*basis, dp)), dp)
    end do
    call assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d_vjp( &
        vertices, triangles, parameters, major_radius, minor_radius, direction, &
        polarization, wave_number, derivative_quadrature_degree, &
        right_hand_side_bar, vertices_bar, parameters_bar, major_radius_bar, &
        minor_radius_bar, direction_bar, polarization_bar, wave_number_bar, status)
    lhs = real(sum(conjg(right_hand_side_bar)*right_hand_side_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + sum(parameters_bar*parameters_dot) + &
        major_radius_bar*major_radius_dot + minor_radius_bar*minor_radius_dot + &
        dot_product(direction_bar, direction_dot) + &
        real(sum(conjg(polarization_bar)*polarization_dot), dp) + &
        wave_number_bar*wave_number_dot
    adjoint_error = abs(lhs - rhs)
    call record_condition(status == 0 .and. adjoint_error < &
        5.0e-9_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "exact-torus plane-wave RHS geometry/data VJP satisfies the adjoint identity")

    allocate(coefficients_dot(size(coefficients)))
    do basis = 1, size(coefficients_dot)
        coefficients_dot(basis) = cmplx( &
            0.013_dp*sin(real(2*basis, dp)), &
            -0.017_dp*cos(real(3*basis, dp)), dp)
    end do
    impedance_dot = -0.023_dp
    call evaluate_maxwell_torus_curved_far_field_rwg_3d_jvp( &
        vertices, triangles, parameters, major_radius, minor_radius, coefficients, &
        direction, wave_number, impedance, derivative_quadrature_degree, &
        vertices_dot, parameters_dot, major_radius_dot, minor_radius_dot, &
        coefficients_dot, direction_dot, wave_number_dot, impedance_dot, far_field, &
        far_field_dot, status)
    call evaluate_maxwell_torus_curved_far_field_rwg_3d( &
        vertices + step*vertices_dot, triangles, parameters + step*parameters_dot, &
        major_radius + step*major_radius_dot, minor_radius + step*minor_radius_dot, &
        coefficients + step*coefficients_dot, &
        rotate_vector(direction, direction_axis, step), wave_number + step*wave_number_dot, &
        impedance + step*impedance_dot, derivative_quadrature_degree, far_field_plus, &
        status)
    call evaluate_maxwell_torus_curved_far_field_rwg_3d( &
        vertices - step*vertices_dot, triangles, parameters - step*parameters_dot, &
        major_radius - step*major_radius_dot, minor_radius - step*minor_radius_dot, &
        coefficients - step*coefficients_dot, &
        rotate_vector(direction, direction_axis, -step), wave_number - step*wave_number_dot, &
        impedance - step*impedance_dot, derivative_quadrature_degree, far_field_minus, &
        status)
    jvp_error = maxval(abs(far_field_dot - &
        (far_field_plus - far_field_minus)/(2.0_dp*step)))
    call record_condition(status == 0 .and. jvp_error < 2.0e-7_dp, &
        "exact-torus Maxwell far-field geometry/data JVP matches reassembly")

    far_field_bar = [ &
        cmplx(0.37_dp, -0.21_dp, dp), cmplx(-0.14_dp, 0.29_dp, dp), &
        cmplx(0.23_dp, 0.11_dp, dp)]
    call evaluate_maxwell_torus_curved_far_field_rwg_3d_vjp( &
        vertices, triangles, parameters, major_radius, minor_radius, coefficients, &
        direction, wave_number, impedance, derivative_quadrature_degree, &
        far_field_bar, vertices_bar, parameters_bar, major_radius_bar, &
        minor_radius_bar, coefficients_bar, direction_bar, wave_number_bar, &
        impedance_bar, status)
    lhs = real(sum(conjg(far_field_bar)*far_field_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + sum(parameters_bar*parameters_dot) + &
        major_radius_bar*major_radius_dot + minor_radius_bar*minor_radius_dot + &
        real(sum(conjg(coefficients_bar)*coefficients_dot), dp) + &
        dot_product(direction_bar, direction_dot) + wave_number_bar*wave_number_dot + &
        impedance_bar*impedance_dot
    adjoint_error = abs(lhs - rhs)
    call record_condition(status == 0 .and. adjoint_error < &
        5.0e-9_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "exact-torus Maxwell far-field geometry/data VJP satisfies the adjoint identity")

    call check_summary("Exact-curved torus Maxwell wave traces")
    if (.not. all_passed) error stop 1

contains

    pure function real_cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function real_cross_product

    pure function rotate_vector(vector, axis, angle) result(rotated)
        real(dp), intent(in) :: vector(3), axis(3), angle
        real(dp) :: rotated(3)

        rotated = cos(angle)*vector + sin(angle)*real_cross_product(axis, vector) + &
            (1.0_dp - cos(angle))*axis*dot_product(axis, vector)
    end function rotate_vector

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_wave
