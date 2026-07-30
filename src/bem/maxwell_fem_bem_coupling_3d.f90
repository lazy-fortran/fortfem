module fortfem_maxwell_fem_bem_coupling_3d
    !! Maxwell trace duality and mixed coupling between the exterior RWG EFIE
    !! operator and arbitrary-order tetrahedral Nedelec fields.
    use fortfem_kinds, only: dp
    use fortfem_assembly_tetra_nedelec_3d, only: &
        assemble_tetra_nedelec_curl_mass_element
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, tetra_nedelec_dof_count, &
        tetra_nedelec_first_kind_t
    use fortfem_tetra_nedelec_global_dof_map, only: &
        build_tetra_nedelec_basis_transform, build_tetra_nedelec_dof_map
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_maxwell_efie_rwg_3d, only: assemble_maxwell_efie_rwg_3d
    use fortfem_maxwell_rwg_surface, only: &
        assemble_maxwell_rwg_mass_matrix, &
        build_maxwell_rwg_surface_space, &
        evaluate_maxwell_rwg_basis, &
        map_maxwell_rwg_to_tetra_nedelec_edges
    use fortfem_tetra_edge_dof_map, only: build_tetra_edge_dof_map
    use fortnum_linalg, only: dense_solve, inv3
    implicit none
    private

    public :: assemble_maxwell_fem_bem_boundary_matrix_3d
    public :: assemble_maxwell_rwg_nedelec_coupling_3d
    public :: assemble_maxwell_fem_bem_system_3d
    public :: solve_maxwell_fem_bem_system_3d

contains

    subroutine assemble_maxwell_rwg_nedelec_coupling_3d( &
            vertices, tetrahedra, boundary_triangles, order, &
            quadrature_degree, coupling, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), boundary_triangles(:, :)
        integer, intent(in) :: order, quadrature_degree
        real(dp), allocatable, intent(out) :: coupling(:, :)
        integer, intent(out) :: status

        type(tetra_nedelec_first_kind_t) :: basis
        integer, allocatable :: edge_orientations(:, :), edges(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :), rwg_edge_triangles(:, :)
        integer, allocatable :: rwg_edges(:, :)
        real(dp), allocatable :: basis_transform(:, :), eta(:), xi(:)
        real(dp), allocatable :: oriented_values(:, :), physical_curls(:, :)
        real(dp), allocatable :: physical_values(:, :), reference_curls(:, :)
        real(dp), allocatable :: reference_values(:, :), weights(:)
        real(dp) :: divergence, inverse(3, 3), jacobian(3, 3), normal(3)
        real(dp) :: panel_jacobian, point(3), reference(3), rwg_value(3)
        integer :: dof, dof_count, inverse_status, node, panel, row, tetrahedron

        status = 1
        if (size(vertices, 1) /= 3 .or. size(tetrahedra, 1) /= 4) return
        if (size(boundary_triangles, 1) /= 3 .or. order < 1) return
        if (quadrature_degree < 0) return
        call initialize_tetra_nedelec_first_kind(order, basis, status)
        if (status /= 0) return
        call build_tetra_nedelec_dof_map( &
            order, tetrahedra, edges, faces, global_dofs, edge_orientations, &
            face_permutations, status)
        if (status /= 0) return
        call build_maxwell_rwg_surface_space( &
            vertices, boundary_triangles, rwg_edges, rwg_edge_triangles, status)
        if (status /= 0 .or. size(rwg_edges, 2) == 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        dof_count = tetra_nedelec_dof_count(basis)
        allocate( &
            coupling(size(rwg_edges, 2), maxval(global_dofs)), &
            basis_transform(dof_count, dof_count), &
            reference_values(3, dof_count), reference_curls(3, dof_count), &
            physical_values(3, dof_count), physical_curls(3, dof_count), &
            oriented_values(3, dof_count))
        coupling = 0.0_dp
        do panel = 1, size(boundary_triangles, 2)
            tetrahedron = adjacent_tetrahedron( &
                tetrahedra, boundary_triangles(:, panel))
            if (tetrahedron == 0) return
            jacobian(:, 1) = vertices(:, tetrahedra(2, tetrahedron)) - &
                vertices(:, tetrahedra(1, tetrahedron))
            jacobian(:, 2) = vertices(:, tetrahedra(3, tetrahedron)) - &
                vertices(:, tetrahedra(1, tetrahedron))
            jacobian(:, 3) = vertices(:, tetrahedra(4, tetrahedron)) - &
                vertices(:, tetrahedra(1, tetrahedron))
            call inv3(jacobian, inverse, inverse_status)
            if (inverse_status /= 0) return
            panel_jacobian = norm2(vector_cross_product( &
                vertices(:, boundary_triangles(2, panel)) - &
                vertices(:, boundary_triangles(1, panel)), &
                vertices(:, boundary_triangles(3, panel)) - &
                vertices(:, boundary_triangles(1, panel))))
            normal = vector_cross_product( &
                vertices(:, boundary_triangles(2, panel)) - &
                vertices(:, boundary_triangles(1, panel)), &
                vertices(:, boundary_triangles(3, panel)) - &
                vertices(:, boundary_triangles(1, panel)))/panel_jacobian
            call build_tetra_nedelec_basis_transform( &
                order, edge_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, status)
            if (status /= 0) return
            do node = 1, size(weights)
                point = &
                    (1.0_dp - xi(node) - eta(node))* &
                    vertices(:, boundary_triangles(1, panel)) + &
                    xi(node)*vertices(:, boundary_triangles(2, panel)) + &
                    eta(node)*vertices(:, boundary_triangles(3, panel))
                reference = matmul( &
                    inverse, point - vertices(:, tetrahedra(1, tetrahedron)))
                call evaluate_tetra_nedelec_first_kind( &
                    basis, reference, reference_values, reference_curls, status)
                if (status /= 0) return
                call map_tetra_nedelec_covariant( &
                    jacobian, reference_values, reference_curls, &
                    physical_values, physical_curls, status)
                if (status /= 0) return
                oriented_values = matmul(physical_values, basis_transform)
                do row = 1, size(rwg_edges, 2)
                    if (.not. any(rwg_edge_triangles(:, row) == panel)) cycle
                    call evaluate_maxwell_rwg_basis( &
                        vertices, boundary_triangles, rwg_edges, &
                        rwg_edge_triangles, row, panel, point, rwg_value, &
                        divergence, status)
                    if (status /= 0) return
                    do dof = 1, dof_count
                        coupling(row, global_dofs(dof, tetrahedron)) = &
                            coupling(row, global_dofs(dof, tetrahedron)) + &
                            panel_jacobian*weights(node)*dot_product( &
                            rwg_value, vector_cross_product( &
                            oriented_values(:, dof), normal))
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_rwg_nedelec_coupling_3d

    subroutine solve_maxwell_fem_bem_system_3d( &
            vertices, tetrahedra, boundary_triangles, curl_coefficient, &
            mass_coefficient, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, right_hand_side, field, current, status, &
            order)
        real(dp), intent(in) :: vertices(:, :), curl_coefficient
        real(dp), intent(in) :: mass_coefficient, wave_number, impedance
        integer, intent(in) :: tetrahedra(:, :), boundary_triangles(:, :)
        integer, intent(in) :: quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), intent(in) :: right_hand_side(:)
        complex(dp), allocatable, intent(out) :: field(:), current(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: order

        complex(dp), allocatable :: matrix(:, :), solution(:)
        integer, allocatable :: edge_orientations(:, :), edges(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :)
        integer :: field_count, info
        integer :: requested_order

        status = 1
        requested_order = 1
        if (present(order)) requested_order = order
        call assemble_maxwell_fem_bem_system_3d( &
            vertices, tetrahedra, boundary_triangles, curl_coefficient, &
            mass_coefficient, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, matrix, status, requested_order)
        if (status /= 0 .or. size(right_hand_side) /= size(matrix, 1)) return
        allocate(solution(size(matrix, 1)))
        call dense_solve(matrix, right_hand_side, solution, info)
        if (info /= 0) then
            status = 2
            return
        end if
        call build_tetra_nedelec_dof_map( &
            requested_order, tetrahedra, edges, faces, global_dofs, &
            edge_orientations, face_permutations, status)
        if (status /= 0) return
        field_count = maxval(global_dofs)
        if (field_count < 1 .or. field_count >= size(matrix, 1)) return
        allocate(field(field_count), current(size(matrix, 1) - field_count))
        field = solution(:field_count)
        current = solution(field_count + 1:)
        status = 0
    end subroutine solve_maxwell_fem_bem_system_3d

    subroutine assemble_maxwell_fem_bem_system_3d( &
            vertices, tetrahedra, boundary_triangles, curl_coefficient, &
            mass_coefficient, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, matrix, status, order)
        real(dp), intent(in) :: vertices(:, :), curl_coefficient
        real(dp), intent(in) :: mass_coefficient, wave_number, impedance
        integer, intent(in) :: tetrahedra(:, :), boundary_triangles(:, :)
        integer, intent(in) :: quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status
        integer, intent(in), optional :: order

        complex(dp), allocatable :: efie(:, :)
        integer, allocatable :: edge_orientations(:, :), edges(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :)
        integer, allocatable :: rwg_edge_triangles(:, :), rwg_edges(:, :)
        real(dp), allocatable :: basis_transform(:, :), coupling(:, :)
        real(dp), allocatable :: element_matrix(:, :), oriented_matrix(:, :)
        real(dp) :: element_vertices(3, 4)
        integer :: column, dof_count, field_count, local_column, local_row
        integer :: node, requested_order, row, rwg_count, tetrahedron

        status = 1
        requested_order = 1
        if (present(order)) requested_order = order
        call build_tetra_nedelec_dof_map( &
            requested_order, tetrahedra, edges, faces, global_dofs, &
            edge_orientations, face_permutations, status)
        if (status /= 0) return
        call build_maxwell_rwg_surface_space( &
            vertices, boundary_triangles, rwg_edges, rwg_edge_triangles, status)
        if (status /= 0 .or. size(rwg_edges, 2) == 0) return
        call assemble_maxwell_rwg_nedelec_coupling_3d( &
            vertices, tetrahedra, boundary_triangles, requested_order, &
            quadrature_degree, coupling, status)
        if (status /= 0) return
        call assemble_maxwell_efie_rwg_3d( &
            vertices, boundary_triangles, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, efie, status)
        if (status /= 0) return
        field_count = maxval(global_dofs)
        rwg_count = size(rwg_edges, 2)
        dof_count = size(global_dofs, 1)
        allocate(matrix(field_count + rwg_count, field_count + rwg_count))
        allocate( &
            basis_transform(dof_count, dof_count), &
            oriented_matrix(dof_count, dof_count))
        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                element_vertices(:, node) = &
                    vertices(:, tetrahedra(node, tetrahedron))
            end do
            call assemble_tetra_nedelec_curl_mass_element( &
                element_vertices, requested_order, quadrature_degree, &
                element_matrix, status, &
                curl_coefficient=curl_coefficient, &
                mass_coefficient=mass_coefficient)
            if (status /= 0) return
            call build_tetra_nedelec_basis_transform( &
                requested_order, edge_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, status)
            if (status /= 0) return
            oriented_matrix = matmul( &
                transpose(basis_transform), &
                matmul(element_matrix, basis_transform))
            do local_column = 1, dof_count
                column = global_dofs(local_column, tetrahedron)
                do local_row = 1, dof_count
                    row = global_dofs(local_row, tetrahedron)
                    matrix(row, column) = matrix(row, column) + &
                        oriented_matrix(local_row, local_column)
                end do
            end do
        end do
        matrix(:field_count, field_count + 1:) = transpose(coupling)
        matrix(field_count + 1:, :field_count) = coupling
        matrix(field_count + 1:, field_count + 1:) = -efie
        status = 0
    end subroutine assemble_maxwell_fem_bem_system_3d

    subroutine assemble_maxwell_fem_bem_boundary_matrix_3d( &
            vertices, tetrahedra, boundary_triangles, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, matrix, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, impedance
        integer, intent(in) :: tetrahedra(:, :), boundary_triangles(:, :)
        integer, intent(in) :: quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: efie(:, :)
        integer, allocatable :: edge_dofs(:, :), edge_orientations(:, :)
        integer, allocatable :: edges(:, :), rwg_dofs(:)
        integer, allocatable :: rwg_edge_triangles(:, :)
        integer, allocatable :: rwg_edges(:, :)
        real(dp), allocatable :: trace_scales(:)
        integer :: column, row

        status = 1
        call build_tetra_edge_dof_map( &
            tetrahedra, edges, edge_dofs, edge_orientations, status)
        if (status /= 0) return
        call build_maxwell_rwg_surface_space( &
            vertices, boundary_triangles, rwg_edges, rwg_edge_triangles, status)
        if (status /= 0 .or. size(rwg_edges, 2) == 0) return
        call map_maxwell_rwg_to_tetra_nedelec_edges( &
            vertices, tetrahedra, rwg_edges, rwg_dofs, trace_scales, status)
        if (status /= 0) return
        call assemble_maxwell_efie_rwg_3d( &
            vertices, boundary_triangles, wave_number, impedance, &
            quadrature_degree, tolerance, max_depth, efie, status)
        if (status /= 0) return
        allocate(matrix(size(edges, 2), size(edges, 2)))
        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        do column = 1, size(rwg_dofs)
            do row = 1, size(rwg_dofs)
                matrix(rwg_dofs(row), rwg_dofs(column)) = &
                    matrix(rwg_dofs(row), rwg_dofs(column)) + &
                    trace_scales(row)*trace_scales(column)*efie(row, column)
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_fem_bem_boundary_matrix_3d

    pure integer function adjacent_tetrahedron( &
            tetrahedra, triangle) result(index)
        integer, intent(in) :: tetrahedra(:, :), triangle(3)

        integer :: tetrahedron, vertex
        logical :: contains

        index = 0
        do tetrahedron = 1, size(tetrahedra, 2)
        contains = .true.
            do vertex = 1, 3
            contains = contains .and. &
                    any(tetrahedra(:, tetrahedron) == triangle(vertex))
            end do
            if (.not. contains) cycle
            index = tetrahedron
            return
        end do
    end function adjacent_tetrahedron

    pure function vector_cross_product(first, second) result(value)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: value(3)

        value = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function vector_cross_product

end module fortfem_maxwell_fem_bem_coupling_3d
