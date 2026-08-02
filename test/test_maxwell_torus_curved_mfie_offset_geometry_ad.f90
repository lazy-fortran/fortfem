program test_maxwell_torus_curved_mfie_offset_geometry_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d, &
        assemble_maxwell_torus_mfie_offset_geometry_jvp, &
        assemble_maxwell_torus_mfie_offset_geometry_vjp, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.4_dp, relative_offset = 0.12_dp
    real(dp), parameter :: step = 2.0e-6_dp
    real(dp) :: major_radius_dot, minor_radius_dot, wave_number_dot
    real(dp) :: relative_offset_dot, major_radius_bar, minor_radius_bar
    real(dp) :: wave_number_bar, relative_offset_bar
    real(dp) :: lhs, rhs, jvp_error
    complex(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    complex(dp), allocatable :: matrix_minus(:, :), matrix_plus(:, :)
    complex(dp), allocatable :: matrix_bar(:, :)
    real(dp), allocatable :: vertices(:, :), parameters(:, :)
    real(dp), allocatable :: vertices_dot(:, :), parameters_dot(:, :)
    real(dp), allocatable :: vertices_minus(:, :), vertices_plus(:, :)
    real(dp), allocatable :: parameters_minus(:, :), parameters_plus(:, :)
    real(dp), allocatable :: vertices_direction(:, :), parameters_direction(:, :)
    real(dp), allocatable :: vertices_bar(:, :), parameters_bar(:, :)
    integer, allocatable :: triangles(:, :)
    integer :: i, j, status, status_minus, status_plus, vertex
    logical :: all_passed

    all_passed = .true.
    major_radius_dot = 0.017_dp
    minor_radius_dot = -0.009_dp
    wave_number_dot = 0.023_dp
    relative_offset_dot = -0.037_dp
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 3, vertices, triangles, parameters)
    allocate( &
        vertices_dot(size(vertices, 1), size(vertices, 2)), &
        parameters_dot(size(parameters, 1), size(parameters, 2)), &
        vertices_minus(size(vertices, 1), size(vertices, 2)), &
        vertices_plus(size(vertices, 1), size(vertices, 2)), &
        parameters_minus(size(parameters, 1), size(parameters, 2)), &
        parameters_plus(size(parameters, 1), size(parameters, 2)), &
        vertices_direction(size(vertices, 1), size(vertices, 2)), &
        parameters_direction(size(parameters, 1), size(parameters, 2)))
    do vertex = 1, size(vertices, 2)
        vertices_dot(:, vertex) = [ &
            0.002_dp*sin(real(vertex, dp)), &
            -0.001_dp*cos(real(2*vertex, dp)), &
            0.0015_dp*sin(real(3*vertex, dp))]
        parameters_dot(:, vertex) = [ &
            0.002_dp*cos(real(vertex, dp)), &
            -0.001_dp*sin(real(2*vertex, dp))]
    end do
    vertices_minus = vertices - step*vertices_dot
    vertices_plus = vertices + step*vertices_dot
    parameters_minus = parameters - step*parameters_dot
    parameters_plus = parameters + step*parameters_dot
    vertices_direction = vertices_dot
    parameters_direction = parameters_dot

    call assemble_maxwell_torus_mfie_offset_geometry_jvp( &
        vertices, triangles, parameters, major_radius, minor_radius, wave_number, &
        2, relative_offset, vertices_dot, parameters_dot, major_radius_dot, &
        minor_radius_dot, wave_number_dot, relative_offset_dot, matrix, &
        matrix_dot, status)
    call assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d( &
        vertices_plus, triangles, parameters_plus, &
        major_radius + step*major_radius_dot, minor_radius + step*minor_radius_dot, &
        wave_number + step*wave_number_dot, 2, &
        relative_offset + step*relative_offset_dot, matrix_plus, status_plus)
    call assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d( &
        vertices_minus, triangles, parameters_minus, &
        major_radius - step*major_radius_dot, minor_radius - step*minor_radius_dot, &
        wave_number - step*wave_number_dot, 2, &
        relative_offset - step*relative_offset_dot, matrix_minus, status_minus)
    jvp_error = maxval(abs(matrix_dot - &
        (matrix_plus - matrix_minus)/(2.0_dp*step)))
    call record_condition(status == 0, &
        "toroidal MFIE geometry JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed toroidal MFIE geometry assembles")
    call record_condition(jvp_error < 2.0e-6_dp, &
        "toroidal MFIE geometry JVP matches reassembly")

    allocate( &
        matrix_bar(size(matrix, 1), size(matrix, 2)), &
        vertices_bar(size(vertices, 1), size(vertices, 2)), &
        parameters_bar(size(parameters, 1), size(parameters, 2)))
    do j = 1, size(matrix_bar, 2)
        do i = 1, size(matrix_bar, 1)
            matrix_bar(i, j) = cmplx( &
                sin(real(i + 2*j, dp)), cos(real(2*i - j, dp)), dp)
        end do
    end do
    call assemble_maxwell_torus_mfie_offset_geometry_vjp( &
        vertices, triangles, parameters, major_radius, minor_radius, wave_number, &
        2, relative_offset, matrix_bar, vertices_bar, parameters_bar, &
        major_radius_bar, minor_radius_bar, wave_number_bar, relative_offset_bar, &
        status)
    lhs = real(sum(conjg(matrix_bar)*matrix_dot), dp)
    rhs = sum(vertices_bar*vertices_direction) + &
        sum(parameters_bar*parameters_direction) + &
        major_radius_bar*major_radius_dot + minor_radius_bar*minor_radius_dot + &
        wave_number_bar*wave_number_dot + relative_offset_bar*relative_offset_dot
    call record_condition(status == 0, &
        "toroidal MFIE geometry VJP succeeds")
    call record_condition(abs(lhs - rhs) < &
        2.0e-7_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "toroidal MFIE geometry products obey the adjoint identity")

    call check_summary("Differentiable toroidal MFIE geometry")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_mfie_offset_geometry_ad
