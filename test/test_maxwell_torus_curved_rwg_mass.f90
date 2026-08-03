program test_maxwell_torus_curved_rwg_mass
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_maxwell_torus_curved_rwg_mass_matrix, &
        assemble_maxwell_torus_curved_rwg_mass_matrix_jvp, &
        assemble_maxwell_torus_curved_rwg_mass_matrix_vjp
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_feec, only: &
        build_maxwell_rwg_surface_space, evaluate_maxwell_torus_curved_rwg_basis, &
        triangle_duffy_quadrature
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: eta(:), matrix(:, :), parameters(:, :)
    real(dp), allocatable :: scaled_matrix(:, :), scaled_vertices(:, :)
    real(dp), allocatable :: values(:), vertices(:, :), weights(:), xi(:)
    real(dp), allocatable :: matrix_dot(:, :), matrix_bar(:, :)
    real(dp), allocatable :: matrix_plus(:, :), matrix_minus(:, :)
    real(dp), allocatable :: vertices_dot(:, :), parameters_dot(:, :)
    real(dp), allocatable :: vertices_bar(:, :), parameters_bar(:, :)
    real(dp) :: divergence, divergence_integral, energy, jacobian
    real(dp) :: normal(3), point(3), relative_scaling_error, tangency_error
    real(dp) :: value(3)
    real(dp) :: major_radius_dot, minor_radius_dot, major_radius_bar
    real(dp) :: minor_radius_bar, step, jvp_error, lhs, rhs, adjoint_error
    integer, parameter :: derivative_quadrature_degree = 4
    integer :: basis, node, panel, status
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 5, 6, vertices, triangles, parameters)
    call build_maxwell_rwg_surface_space( &
        vertices, triangles, edge_vertices, edge_triangles, status)
    call record_condition(status == 0, "toroidal RWG topology builds")
    call triangle_duffy_quadrature(12, xi, eta, weights, status)
    basis = 1
    divergence_integral = 0.0_dp
    tangency_error = 0.0_dp
    do panel = 1, size(triangles, 2)
        if (.not. any(edge_triangles(:, basis) == panel)) cycle
        do node = 1, size(weights)
            call evaluate_maxwell_torus_curved_rwg_basis( &
                vertices, triangles, parameters, edge_vertices, &
                edge_triangles, basis, panel, major_radius, minor_radius, &
                xi(node), eta(node), point, value, divergence, jacobian, status)
            normal = torus_normal(point)
            tangency_error = max( &
                tangency_error, abs(dot_product(normal, value)))
            divergence_integral = divergence_integral + &
                weights(node)*jacobian*divergence
        end do
    end do
    call record_condition(status == 0 .and. tangency_error < 3.0e-14_dp, &
        "surface-Piola RWG values are tangent to the exact torus")
    call record_condition(abs(divergence_integral) < 3.0e-14_dp, &
        "interior RWG divergence has zero closed-surface integral")

    call assemble_maxwell_torus_curved_rwg_mass_matrix( &
        vertices, triangles, parameters, major_radius, minor_radius, 12, &
        matrix, status)
    allocate(values(size(matrix, 1)))
    do basis = 1, size(values)
        values(basis) = sin(real(3*basis + 1, dp))
    end do
    energy = dot_product(values, matmul(matrix, values))
    call record_condition(status == 0 .and. &
        maxval(abs(matrix - transpose(matrix))) < 3.0e-14_dp, &
        "curved-torus RWG mass is symmetric")
    call record_condition(energy > 0.0_dp, &
        "curved-torus RWG mass has positive field energy")

    allocate( &
        vertices_dot(size(vertices, 1), size(vertices, 2)), &
        parameters_dot(size(parameters, 1), size(parameters, 2)))
    do basis = 1, size(vertices_dot, 2)
        vertices_dot(:, basis) = [ &
            0.013_dp*sin(real(2*basis, dp)), &
            -0.009_dp*cos(real(3*basis, dp)), &
            0.011_dp*sin(real(5*basis, dp))]
        parameters_dot(:, basis) = [ &
            0.007_dp*cos(real(basis + 1, dp)), &
            -0.005_dp*sin(real(2*basis + 1, dp))]
    end do
    major_radius_dot = 0.021_dp
    minor_radius_dot = -0.013_dp
    call assemble_maxwell_torus_curved_rwg_mass_matrix_jvp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        derivative_quadrature_degree, &
        vertices_dot, parameters_dot, major_radius_dot, minor_radius_dot, &
        matrix, matrix_dot, status)
    step = 2.0e-6_dp
    call assemble_maxwell_torus_curved_rwg_mass_matrix( &
        vertices + step*vertices_dot, triangles, parameters + step*parameters_dot, &
        major_radius + step*major_radius_dot, minor_radius + step*minor_radius_dot, &
        derivative_quadrature_degree, matrix_plus, status)
    call assemble_maxwell_torus_curved_rwg_mass_matrix( &
        vertices - step*vertices_dot, triangles, parameters - step*parameters_dot, &
        major_radius - step*major_radius_dot, minor_radius - step*minor_radius_dot, &
        derivative_quadrature_degree, matrix_minus, status)
    jvp_error = maxval(abs(matrix_dot - (matrix_plus - matrix_minus)/(2.0_dp*step)))
    call record_condition(status == 0 .and. jvp_error < 2.0e-7_dp, &
        "curved-torus RWG mass geometry JVP matches reassembly")

    allocate(matrix_bar(size(matrix, 1), size(matrix, 2)))
    do basis = 1, size(matrix_bar, 2)
        do panel = 1, size(matrix_bar, 1)
            matrix_bar(panel, basis) = sin(real(2*panel + 3*basis, dp))
        end do
    end do
    allocate( &
        vertices_bar(size(vertices, 1), size(vertices, 2)), &
        parameters_bar(size(parameters, 1), size(parameters, 2)))
    call assemble_maxwell_torus_curved_rwg_mass_matrix_vjp( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        derivative_quadrature_degree, matrix_bar, &
        vertices_bar, parameters_bar, major_radius_bar, minor_radius_bar, status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(vertices_bar*vertices_dot) + sum(parameters_bar*parameters_dot) + &
        major_radius_bar*major_radius_dot + minor_radius_bar*minor_radius_dot
    adjoint_error = abs(lhs - rhs)
    call record_condition(status == 0 .and. adjoint_error < &
        5.0e-9_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "curved-torus RWG mass geometry VJP satisfies the adjoint identity")

    scaled_vertices = 3.0_dp*vertices
    call assemble_maxwell_torus_curved_rwg_mass_matrix( &
        scaled_vertices, triangles, parameters, 3.0_dp*major_radius, &
        3.0_dp*minor_radius, 12, scaled_matrix, status)
    call assemble_maxwell_torus_curved_rwg_mass_matrix( &
        vertices, triangles, parameters, major_radius, minor_radius, 12, &
        matrix, status)
    relative_scaling_error = maxval(abs(scaled_matrix - 9.0_dp*matrix))/ &
        maxval(abs(scaled_matrix))
    call record_condition(status == 0 .and. relative_scaling_error < &
        5.0e-14_dp, "curved-torus RWG mass obeys analytical area scaling")

    call check_summary("Exact-curved torus RWG mass matrix")
    if (.not. all_passed) error stop 1

contains

    pure function torus_normal(position) result(normal_value)
        real(dp), intent(in) :: position(3)
        real(dp) :: normal_value(3)
        real(dp) :: cylindrical_radius

        cylindrical_radius = sqrt(position(1)**2 + position(2)**2)
        normal_value = [ &
            (cylindrical_radius - major_radius)*position(1)/ &
            cylindrical_radius, &
            (cylindrical_radius - major_radius)*position(2)/ &
            cylindrical_radius, position(3)]/minor_radius
    end function torus_normal

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_rwg_mass
