program test_maxwell_sphere_curved_magnetic_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d, &
        evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d_jvp, &
        evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d_vjp, &
        generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: radius = 1.0_dp
    real(dp), parameter :: observation(3) = [2.2_dp, 0.3_dp, -0.1_dp]
    real(dp), parameter :: wave_number = 0.7_dp
    real(dp), parameter :: step = 1.0e-6_dp
    real(dp), allocatable :: vertices(:, :), vertices_dot(:, :)
    real(dp), allocatable :: vertices_bar(:, :)
    integer, allocatable :: triangles(:, :)
    complex(dp), allocatable :: coefficients(:), coefficients_dot(:)
    complex(dp), allocatable :: coefficients_bar(:)
    complex(dp) :: field_dot(3), field_bar(3), field_plus(3), field_minus(3)
    real(dp) :: radius_dot, radius_bar, wave_number_dot, wave_number_bar
    real(dp) :: observation_dot(3), observation_bar(3)
    real(dp) :: derivative_error, adjoint_error, lhs, rhs
    integer :: i, status
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(radius, 0, vertices, triangles)
    allocate(vertices_dot(size(vertices, 1), size(vertices, 2)))
    vertices_dot = 0.0_dp
    do i = 1, size(vertices, 2)
        vertices_dot(:, i) = [ &
            0.002_dp*i, -0.0013_dp*i, 0.0008_dp*i]
    end do
    call build_coefficients(coefficients, coefficients_dot, status)
    call record_condition(status == 0, &
        "curved-sphere magnetic derivative coefficients initialize")
    radius_dot = 0.003_dp
    observation_dot = [0.004_dp, -0.003_dp, 0.002_dp]
    wave_number_dot = 0.017_dp

    call evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d_jvp( &
        vertices, triangles, radius, coefficients, observation, wave_number, 8, &
        vertices_dot, radius_dot, coefficients_dot, observation_dot, &
        wave_number_dot, field_dot, status)
    call evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d( &
        vertices + step*vertices_dot, triangles, radius + step*radius_dot, &
        coefficients + step*coefficients_dot, observation + step*observation_dot, &
        wave_number + step*wave_number_dot, 8, field_plus, status)
    call evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d( &
        vertices - step*vertices_dot, triangles, radius - step*radius_dot, &
        coefficients - step*coefficients_dot, observation - step*observation_dot, &
        wave_number - step*wave_number_dot, 8, field_minus, status)
    derivative_error = maxval(abs(field_dot - &
        (field_plus - field_minus)/(2.0_dp*step)))
    call record_condition(status == 0 .and. derivative_error < 2.0e-7_dp, &
        "curved-sphere magnetic JVP matches central reassembly")

    field_bar = [ &
        cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.1_dp, 0.6_dp, dp)]
    allocate(vertices_bar(size(vertices, 1), size(vertices, 2)))
    call evaluate_maxwell_sphere_curved_magnetic_field_rwg_3d_vjp( &
        vertices, triangles, radius, coefficients, observation, wave_number, 8, &
        field_bar, vertices_bar, radius_bar, coefficients_bar, observation_bar, &
        wave_number_bar, status)
    lhs = real(sum(conjg(field_bar)*field_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + radius_bar*radius_dot + &
        real(sum(conjg(coefficients_bar)*coefficients_dot), dp) + &
        dot_product(observation_bar, observation_dot) + &
        wave_number_bar*wave_number_dot
    adjoint_error = abs(lhs - rhs)
    call record_condition(status == 0 .and. adjoint_error < &
        5.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "curved-sphere magnetic VJP satisfies the real complex adjoint")

    call check_summary("curved-sphere magnetic derivative")
    if (.not. all_passed) error stop 1

contains

    subroutine build_coefficients(values, values_dot, local_status)
        complex(dp), allocatable, intent(out) :: values(:), values_dot(:)
        integer, intent(out) :: local_status

        integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
        integer :: basis

        call build_maxwell_rwg_surface_space_wrapper( &
            vertices, triangles, edge_vertices, edge_triangles, local_status)
        if (local_status /= 0) return
        allocate(values(size(edge_vertices, 2)), values_dot(size(edge_vertices, 2)))
        do basis = 1, size(values)
            values(basis) = cmplx( &
                cos(0.31_dp*basis), sin(0.17_dp*basis), dp)
            values_dot(basis) = cmplx( &
                0.11_dp*sin(0.13_dp*basis), -0.07_dp*cos(0.19_dp*basis), dp)
        end do
        local_status = 0
    end subroutine build_coefficients

    subroutine build_maxwell_rwg_surface_space_wrapper( &
            surface_vertices, surface_triangles, edge_vertices, edge_triangles, &
            local_status)
        use fortfem_api, only: build_maxwell_rwg_surface_space
        real(dp), intent(in) :: surface_vertices(:, :)
        integer, intent(in) :: surface_triangles(:, :)
        integer, allocatable, intent(out) :: edge_vertices(:, :)
        integer, allocatable, intent(out) :: edge_triangles(:, :)
        integer, intent(out) :: local_status

        call build_maxwell_rwg_surface_space( &
            surface_vertices, surface_triangles, edge_vertices, edge_triangles, &
            local_status)
    end subroutine build_maxwell_rwg_surface_space_wrapper

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_magnetic_ad
