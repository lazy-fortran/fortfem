program test_nedelec_weighted_assembly
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_assembly_nedelec_2d, only: &
        assemble_nedelec_weighted_element, assemble_nedelec_weighted, &
        assemble_nedelec_axisymmetric_fourier
    use check, only: check_condition, check_summary
    implicit none

    logical :: all_passed

    all_passed = .true.
    call test_radial_curl_weight()
    call test_global_radial_mass_energy()
    call test_axisymmetric_fourier_mass()
    call check_summary("Weighted Nedelec assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine test_radial_curl_weight()
        type(mesh_2d_t) :: mesh
        real(dp) :: element_matrix(3, 3)

        mesh%n_vertices = 3
        mesh%n_triangles = 1
        allocate(mesh%vertices(2, 3), mesh%triangles(3, 1))
        mesh%vertices(:, 1) = [0.0_dp, 0.0_dp]
        mesh%vertices(:, 2) = [1.0_dp, 0.0_dp]
        mesh%vertices(:, 3) = [0.0_dp, 1.0_dp]
        mesh%triangles(:, 1) = [1, 2, 3]

        call assemble_nedelec_weighted_element(mesh, 1, &
            radial_coordinate, zero_coefficient, 3, element_matrix)

        call record_condition( &
            maxval(abs(element_matrix - 2.0_dp / 3.0_dp)) < 1.0e-13_dp, &
            "R-weighted curl matrix matches exact linear integration")
    end subroutine test_radial_curl_weight

    subroutine test_global_radial_mass_energy()
        type(mesh_2d_t) :: mesh
        real(dp), allocatable :: matrix(:, :), x_dofs(:), y_dofs(:)
        real(dp), allocatable :: matrix_times_x(:), matrix_times_y(:)
        real(dp) :: x_energy, y_energy, cross_energy
        integer :: edge_id, dof_id, vertex_a, vertex_b

        mesh%n_vertices = 4
        mesh%n_triangles = 2
        mesh%has_triangles = .true.
        allocate(mesh%vertices(2, 4), mesh%triangles(3, 2))
        mesh%vertices(:, 1) = [0.0_dp, 0.0_dp]
        mesh%vertices(:, 2) = [1.0_dp, 0.0_dp]
        mesh%vertices(:, 3) = [1.0_dp, 1.0_dp]
        mesh%vertices(:, 4) = [0.0_dp, 1.0_dp]
        mesh%triangles(:, 1) = [1, 2, 3]
        mesh%triangles(:, 2) = [1, 3, 4]
        call mesh%build_edge_connectivity()

        allocate(matrix(mesh%n_edges, mesh%n_edges))
        allocate(x_dofs(mesh%n_edges), y_dofs(mesh%n_edges))
        allocate(matrix_times_x(mesh%n_edges), matrix_times_y(mesh%n_edges))
        call assemble_nedelec_weighted(mesh, zero_coefficient, &
            radial_coordinate, 3, matrix)

        do edge_id = 1, mesh%n_edges
            vertex_a = mesh%edges(1, edge_id)
            vertex_b = mesh%edges(2, edge_id)
            dof_id = mesh%edge_to_dof(edge_id) + 1
            x_dofs(dof_id) = mesh%vertices(1, vertex_b) - &
                mesh%vertices(1, vertex_a)
            y_dofs(dof_id) = mesh%vertices(2, vertex_b) - &
                mesh%vertices(2, vertex_a)
        end do

        matrix_times_x = matmul(matrix, x_dofs)
        matrix_times_y = matmul(matrix, y_dofs)
        x_energy = dot_product(x_dofs, matrix_times_x)
        y_energy = dot_product(y_dofs, matrix_times_y)
        cross_energy = dot_product(x_dofs, matrix_times_y)

        call record_condition(abs(x_energy - 0.5_dp) < 1.0e-13_dp, &
            "Constant x field has exact integral of R over the square")
        call record_condition(abs(y_energy - 0.5_dp) < 1.0e-13_dp, &
            "Constant y field has exact integral of R over the square")
        call record_condition(abs(cross_energy) < 1.0e-13_dp, &
            "R-weighted orthogonal constant fields have zero cross energy")
    end subroutine test_global_radial_mass_energy

    subroutine test_axisymmetric_fourier_mass()
        type(mesh_2d_t) :: mesh
        real(dp), allocatable :: matrix(:, :), x_dofs(:), matrix_times_x(:)
        real(dp) :: energies(3), errors(3), exact_energy
        integer, parameter :: degrees(3) = [3, 5, 7]
        integer :: edge_id, dof_id, vertex_a, vertex_b, degree_id

        mesh%n_vertices = 4
        mesh%n_triangles = 2
        mesh%has_triangles = .true.
        allocate(mesh%vertices(2, 4), mesh%triangles(3, 2))
        mesh%vertices(:, 1) = [1.0_dp, 0.0_dp]
        mesh%vertices(:, 2) = [2.0_dp, 0.0_dp]
        mesh%vertices(:, 3) = [2.0_dp, 1.0_dp]
        mesh%vertices(:, 4) = [1.0_dp, 1.0_dp]
        mesh%triangles(:, 1) = [1, 2, 3]
        mesh%triangles(:, 2) = [1, 3, 4]
        call mesh%build_edge_connectivity()

        allocate(matrix(mesh%n_edges, mesh%n_edges))
        allocate(x_dofs(mesh%n_edges), matrix_times_x(mesh%n_edges))
        call mesh%build_edge_dof_numbering()
        do edge_id = 1, mesh%n_edges
            vertex_a = mesh%edges(1, edge_id)
            vertex_b = mesh%edges(2, edge_id)
            dof_id = mesh%edge_to_dof(edge_id) + 1
            x_dofs(dof_id) = mesh%vertices(1, vertex_b) - &
                mesh%vertices(1, vertex_a)
        end do

        do degree_id = 1, size(degrees)
            call assemble_nedelec_axisymmetric_fourier( &
                mesh, 2, degrees(degree_id), matrix)
            matrix_times_x = matmul(matrix, x_dofs)
            energies(degree_id) = dot_product(x_dofs, matrix_times_x)
        end do
        exact_energy = 4.0_dp * log(2.0_dp)
        errors = abs(energies - exact_energy)

        call record_condition(errors(2) < errors(1), &
            "Axisymmetric Fourier mass converges from degree three to five")
        call record_condition(errors(3) < errors(2), &
            "Axisymmetric Fourier mass converges from degree five to seven")
        call record_condition(errors(3) < 2.0e-6_dp, &
            "Axisymmetric Fourier mass resolves the logarithmic exact energy")
    end subroutine test_axisymmetric_fourier_mass

    pure real(dp) function radial_coordinate(x, y) result(value)
        real(dp), intent(in) :: x, y

        value = x
        associate (unused_y => y)
            if (kind(unused_y) /= dp) error stop
        end associate
    end function radial_coordinate

    pure real(dp) function zero_coefficient(x, y) result(value)
        real(dp), intent(in) :: x, y

        value = 0.0_dp
        associate (unused_x => x, unused_y => y)
            if (kind(unused_x) /= dp .or. kind(unused_y) /= dp) error stop
        end associate
    end function zero_coefficient

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_nedelec_weighted_assembly
