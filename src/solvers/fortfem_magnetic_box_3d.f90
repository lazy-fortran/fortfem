module fortfem_magnetic_box_3d
    use fortfem_assembly_tetra_nedelec_3d, only: &
        assemble_tetra_nedelec_vector_load, &
        assemble_tetra_nedelec_weighted_csc
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_solve_csc
    use fortfem_structured_tetra_box_mesh, only: &
        generate_structured_tetra_box_mesh
    use fortfem_tetra_edge_dof_map, only: build_tetra_edge_dof_map
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, tetra_nedelec_first_kind_t
    use fortfem_tetra_nedelec_first_order, only: &
        evaluate_tetra_nedelec_first_order
    use fortfem_tetra_nedelec_global_dof_map, only: &
        build_tetra_nedelec_basis_transform, build_tetra_nedelec_dof_map
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant
    use fortfem_tetra_nedelec_solver_3d, only: &
        solve_tetra_nedelec_weighted_curl_mass
    use fortnum_linalg, only: det3, inv3
    use fortsparse, only: csc_from_triplet, csc_t, &
        FORTSPARSE_INTERNAL_ERROR, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none

    private

    real(dp), parameter :: gauge_mass = 1.0e-10_dp
    ! Higher-order curl-curl systems need a resolvable penalty on the gradient
    ! nullspace; this remains six orders below the physical operator scale.
    real(dp), parameter :: higher_order_gauge_mass = 1.0e-6_dp
    real(dp), parameter :: box_bounds(3, 2) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [3, 2])

    public :: solve_magnetic_box_3d

contains

    subroutine solve_magnetic_box_3d( &
            n_xy, n_z, az_centre, n_dofs, status, magnetic_point, &
            magnetic_field, order)
        integer, intent(in) :: n_xy, n_z
        real(dp), intent(out) :: az_centre
        integer, intent(out) :: n_dofs
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: magnetic_point(3)
        real(dp), intent(out), optional :: magnetic_field(3)
        integer, intent(in), optional :: order

        type(csc_t) :: matrix
        integer, allocatable :: edges(:, :), global_dofs(:, :)
        integer, allocatable :: orientations(:, :), tetrahedra(:, :)
        real(dp), allocatable :: right_hand_side(:), solution(:)
        real(dp), allocatable :: vertices(:, :)
        integer :: local_status, polynomial_order

        az_centre = 0.0_dp
        n_dofs = 0
        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Magnetic box solve failed")
        if (present(magnetic_point)) then
            if (.not. present(magnetic_field)) return
        end if
        if (present(magnetic_field)) then
            if (.not. present(magnetic_point)) return
        end if
        polynomial_order = 1
        if (present(order)) polynomial_order = order
        if (polynomial_order < 1) return
        call generate_structured_tetra_box_mesh( &
            box_bounds, [n_xy, n_xy, n_z], vertices, tetrahedra, local_status)
        if (local_status /= 0) return
        if (polynomial_order > 1) then
            call solve_tetra_nedelec_weighted_curl_mass( &
                vertices, tetrahedra, polynomial_order, box_reluctivity, &
                axial_source, higher_order_gauge_mass, solution, status, &
                .true.)
            if (status%code /= 0) return
            n_dofs = size(solution)
            call probe_box_solution_order( &
                vertices, tetrahedra, polynomial_order, solution, &
                [0.5_dp, 0.5_dp, 1.5_dp], az_centre, local_status)
            if (local_status /= 0) then
                call status_set( &
                    status, FORTSPARSE_INTERNAL_ERROR, &
                    "Higher-order magnetic box centre probe failed")
                return
            end if
            if (present(magnetic_point)) then
                call probe_box_curl_order( &
                    vertices, tetrahedra, polynomial_order, solution, &
                    magnetic_point, magnetic_field, local_status)
                if (local_status /= 0) then
                    call status_set( &
                        status, FORTSPARSE_INTERNAL_ERROR, &
                        "Higher-order magnetic box curl probe failed")
                    return
                end if
            end if
            call status_set(status, 0, "")
            return
        end if
        call build_tetra_edge_dof_map( &
            tetrahedra, edges, global_dofs, orientations, local_status)
        if (local_status /= 0) return
        n_dofs = size(edges, 2)
        call assemble_tetra_nedelec_weighted_csc( &
            vertices, tetrahedra, box_reluctivity, gauge_mass, matrix, status)
        if (status%code /= 0) return
        call assemble_tetra_nedelec_vector_load( &
            vertices, tetrahedra, axial_source, right_hand_side, status)
        if (status%code /= 0) return
        call solve_box_interior( &
            vertices, edges, matrix, right_hand_side, solution, status)
        if (status%code /= 0) return
        call probe_box_az( &
            vertices, tetrahedra, global_dofs, orientations, solution, &
            az_centre, local_status)
        if (local_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INTERNAL_ERROR, &
                "Magnetic box centre probe failed")
            return
        end if
        if (present(magnetic_point)) then
            call probe_box_curl( &
                vertices, tetrahedra, global_dofs, orientations, solution, &
                magnetic_point, magnetic_field, local_status)
            if (local_status /= 0) then
                call status_set( &
                    status, FORTSPARSE_INTERNAL_ERROR, &
                    "Magnetic box curl probe failed")
                return
            end if
        end if
        call status_set(status, 0, "")
    end subroutine solve_magnetic_box_3d

    subroutine solve_box_interior( &
            vertices, edges, matrix, rhs, solution, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: edges(:, :)
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: rhs(:)
        real(dp), allocatable, intent(out) :: solution(:)
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: interior_matrix
        integer, allocatable :: columns(:), interior_index(:), rows(:)
        real(dp), allocatable :: interior_rhs(:), interior_solution(:)
        real(dp), allocatable :: values(:)
        integer :: column, edge, entry, interior_count, interior_entry
        integer :: row, solve_status

        allocate(interior_index(size(edges, 2)))
        interior_index = 0
        interior_count = 0
        do edge = 1, size(edges, 2)
            if (.not. box_boundary_edge(vertices, edges(:, edge))) then
                interior_count = interior_count + 1
                interior_index(edge) = interior_count
            end if
        end do
        if (interior_count < 1) then
            call status_set( &
                status, FORTSPARSE_INVALID_MATRIX, &
                "Magnetic box has no interior edge degrees of freedom")
            return
        end if
        allocate(rows(matrix%nnz), columns(matrix%nnz), values(matrix%nnz))
        interior_entry = 0
        do column = 1, matrix%ncol
            if (interior_index(column) == 0) cycle
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                if (interior_index(row) == 0) cycle
                interior_entry = interior_entry + 1
                rows(interior_entry) = interior_index(row)
                columns(interior_entry) = interior_index(column)
                values(interior_entry) = matrix%val(entry)
            end do
        end do
        call csc_from_triplet( &
            interior_count, interior_count, rows(:interior_entry), &
            columns(:interior_entry), values(:interior_entry), &
            interior_matrix, status)
        if (status%code /= 0) return
        allocate(interior_rhs(interior_count))
        allocate(interior_solution(interior_count))
        do edge = 1, size(edges, 2)
            if (interior_index(edge) > 0) then
                interior_rhs(interior_index(edge)) = rhs(edge)
            end if
        end do
        call sparse_direct_solve_csc( &
            interior_count, interior_matrix%col_ptr, interior_matrix%row_idx, &
            interior_matrix%val, interior_rhs, interior_solution, solve_status)
        if (solve_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INTERNAL_ERROR, &
                "Magnetic box sparse solve failed")
            return
        end if
        allocate(solution(size(edges, 2)))
        solution = 0.0_dp
        do edge = 1, size(edges, 2)
            if (interior_index(edge) > 0) then
                solution(edge) = interior_solution(interior_index(edge))
            end if
        end do
        call status_set(status, 0, "")
    end subroutine solve_box_interior

    pure logical function box_boundary_edge(vertices, edge_vertices)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: edge_vertices(2)

        real(dp), parameter :: tolerance = 64.0_dp * epsilon(1.0_dp)
        real(dp) :: first(3), second(3)

        first = vertices(:, edge_vertices(1))
        second = vertices(:, edge_vertices(2))
        box_boundary_edge = &
            same_plane(first(1), second(1), 0.0_dp, tolerance) .or. &
            same_plane(first(1), second(1), 1.0_dp, tolerance) .or. &
            same_plane(first(2), second(2), 0.0_dp, tolerance) .or. &
            same_plane(first(2), second(2), 1.0_dp, tolerance) .or. &
            same_plane(first(3), second(3), 1.0_dp, tolerance) .or. &
            same_plane(first(3), second(3), 2.0_dp, tolerance)
    end function box_boundary_edge

    pure logical function same_plane(first, second, plane, tolerance)
        real(dp), intent(in) :: first, second, plane, tolerance

        same_plane = abs(first - plane) <= tolerance .and. &
            abs(second - plane) <= tolerance
    end function same_plane

    subroutine probe_box_solution_order( &
            vertices, tetrahedra, order, solution, point, az_value, status)
        real(dp), intent(in) :: vertices(:, :), solution(:), point(3)
        integer, intent(in) :: tetrahedra(:, :), order
        real(dp), intent(out) :: az_value
        integer, intent(out) :: status

        real(dp) :: curl_value(3), value(3)

        call evaluate_box_solution_order( &
            vertices, tetrahedra, order, solution, point, value, curl_value, &
            status)
        az_value = value(3)
    end subroutine probe_box_solution_order

    subroutine probe_box_curl_order( &
            vertices, tetrahedra, order, solution, point, curl_value, status)
        real(dp), intent(in) :: vertices(:, :), solution(:), point(3)
        integer, intent(in) :: tetrahedra(:, :), order
        real(dp), intent(out) :: curl_value(3)
        integer, intent(out) :: status

        real(dp) :: value(3)

        call evaluate_box_solution_order( &
            vertices, tetrahedra, order, solution, point, value, curl_value, &
            status)
    end subroutine probe_box_curl_order

    subroutine evaluate_box_solution_order( &
            vertices, tetrahedra, order, solution, point, value, curl_value, &
            status)
        real(dp), intent(in) :: vertices(:, :), solution(:), point(3)
        integer, intent(in) :: tetrahedra(:, :), order
        real(dp), intent(out) :: value(3), curl_value(3)
        integer, intent(out) :: status

        type(tetra_nedelec_first_kind_t) :: basis
        integer, allocatable :: edge_orientations(:, :), edges(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :)
        real(dp), allocatable :: basis_transform(:, :), local_dofs(:)
        real(dp), allocatable :: physical_curls(:, :), physical_values(:, :)
        real(dp), allocatable :: reference_curls(:, :)
        real(dp), allocatable :: reference_values(:, :)
        real(dp) :: determinant, inverse_jacobian(3, 3), jacobian(3, 3)
        real(dp) :: reference_point(3), vertex_a(3)
        integer :: containing_count, dof_count, inverse_status, tetrahedron

        value = 0.0_dp
        curl_value = 0.0_dp
        containing_count = 0
        status = 1
        call build_tetra_nedelec_dof_map( &
            order, tetrahedra, edges, faces, global_dofs, edge_orientations, &
            face_permutations, status)
        if (status /= 0) return
        if (size(solution) /= maxval(global_dofs)) return
        call initialize_tetra_nedelec_first_kind(order, basis, status)
        if (status /= 0) return
        dof_count = size(global_dofs, 1)
        allocate ( &
            basis_transform(dof_count, dof_count), local_dofs(dof_count), &
            reference_values(3, dof_count), reference_curls(3, dof_count), &
            physical_values(3, dof_count), physical_curls(3, dof_count))

        do tetrahedron = 1, size(tetrahedra, 2)
            vertex_a = vertices(:, tetrahedra(1, tetrahedron))
            jacobian(:, 1) = &
                vertices(:, tetrahedra(2, tetrahedron)) - vertex_a
            jacobian(:, 2) = &
                vertices(:, tetrahedra(3, tetrahedron)) - vertex_a
            jacobian(:, 3) = &
                vertices(:, tetrahedra(4, tetrahedron)) - vertex_a
            determinant = det3(jacobian)
            if (determinant <= 0.0_dp) cycle
            call inv3(jacobian, inverse_jacobian, inverse_status)
            if (inverse_status /= 0) cycle
            reference_point = matmul(inverse_jacobian, point - vertex_a)
            if (any(reference_point < -1.0e-12_dp)) cycle
            if (sum(reference_point) > 1.0_dp + 1.0e-12_dp) cycle
            call build_tetra_nedelec_basis_transform( &
                order, edge_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, status)
            if (status /= 0) return
            local_dofs = matmul( &
                basis_transform, solution(global_dofs(:, tetrahedron)))
            call evaluate_tetra_nedelec_first_kind( &
                basis, reference_point, reference_values, reference_curls, &
                status)
            if (status /= 0) return
            call map_tetra_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, status)
            if (status /= 0) return
            value = value + matmul(physical_values, local_dofs)
            curl_value = curl_value + matmul(physical_curls, local_dofs)
            containing_count = containing_count + 1
        end do
        if (containing_count < 1) return
        value = value/real(containing_count, dp)
        curl_value = curl_value/real(containing_count, dp)
        status = 0
    end subroutine evaluate_box_solution_order

    subroutine probe_box_az( &
            vertices, tetrahedra, global_dofs, orientations, solution, &
            az_centre, status)
        real(dp), intent(in) :: vertices(:, :), solution(:)
        integer, intent(in) :: tetrahedra(:, :)
        integer, intent(in) :: global_dofs(:, :), orientations(:, :)
        real(dp), intent(out) :: az_centre
        integer, intent(out) :: status

        real(dp) :: determinant, inverse_jacobian(3, 3), jacobian(3, 3)
        real(dp) :: local_dofs(6), physical_curls(3, 6)
        real(dp) :: physical_values(3, 6), point(3), reference_curls(3, 6)
        real(dp) :: reference_point(3), reference_values(3, 6), value(3)
        real(dp) :: vertex_a(3)
        integer :: dof, inverse_status, tetrahedron

        az_centre = 0.0_dp
        status = 1
        point = [0.5_dp, 0.5_dp, 1.5_dp]
        do tetrahedron = 1, size(tetrahedra, 2)
            vertex_a = vertices(:, tetrahedra(1, tetrahedron))
            jacobian(:, 1) = vertices(:, tetrahedra(2, tetrahedron)) - vertex_a
            jacobian(:, 2) = vertices(:, tetrahedra(3, tetrahedron)) - vertex_a
            jacobian(:, 3) = vertices(:, tetrahedra(4, tetrahedron)) - vertex_a
            determinant = det3(jacobian)
            if (determinant <= 0.0_dp) cycle
            call inv3(jacobian, inverse_jacobian, inverse_status)
            if (inverse_status /= 0) cycle
            reference_point = matmul(inverse_jacobian, point - vertex_a)
            if (any(reference_point < -1.0e-12_dp)) cycle
            if (sum(reference_point) > 1.0_dp + 1.0e-12_dp) cycle
            call evaluate_tetra_nedelec_first_order( &
                reference_point, reference_values, reference_curls, status)
            if (status /= 0) return
            call map_tetra_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, status)
            if (status /= 0) return
            do dof = 1, 6
                local_dofs(dof) = real( &
                    orientations(dof, tetrahedron), dp) * &
                    solution(global_dofs(dof, tetrahedron))
            end do
            value = matmul(physical_values, local_dofs)
            az_centre = value(3)
            status = 0
            return
        end do
    end subroutine probe_box_az

    subroutine probe_box_curl( &
            vertices, tetrahedra, global_dofs, orientations, solution, point, &
            curl_value, status)
        real(dp), intent(in) :: vertices(:, :), solution(:), point(3)
        integer, intent(in) :: tetrahedra(:, :)
        integer, intent(in) :: global_dofs(:, :), orientations(:, :)
        real(dp), intent(out) :: curl_value(3)
        integer, intent(out) :: status

        real(dp) :: determinant, inverse_jacobian(3, 3), jacobian(3, 3)
        real(dp) :: local_dofs(6), physical_curls(3, 6)
        real(dp) :: physical_values(3, 6), reference_curls(3, 6)
        real(dp) :: reference_point(3), reference_values(3, 6), vertex_a(3)
        integer :: dof, inverse_status, tetrahedron

        curl_value = 0.0_dp
        status = 1
        do tetrahedron = 1, size(tetrahedra, 2)
            vertex_a = vertices(:, tetrahedra(1, tetrahedron))
            jacobian(:, 1) = vertices(:, tetrahedra(2, tetrahedron)) - vertex_a
            jacobian(:, 2) = vertices(:, tetrahedra(3, tetrahedron)) - vertex_a
            jacobian(:, 3) = vertices(:, tetrahedra(4, tetrahedron)) - vertex_a
            determinant = det3(jacobian)
            if (determinant <= 0.0_dp) cycle
            call inv3(jacobian, inverse_jacobian, inverse_status)
            if (inverse_status /= 0) cycle
            reference_point = matmul(inverse_jacobian, point - vertex_a)
            if (any(reference_point < -1.0e-12_dp)) cycle
            if (sum(reference_point) > 1.0_dp + 1.0e-12_dp) cycle
            call evaluate_tetra_nedelec_first_order( &
                reference_point, reference_values, reference_curls, status)
            if (status /= 0) return
            call map_tetra_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, status)
            if (status /= 0) return
            do dof = 1, 6
                local_dofs(dof) = real( &
                    orientations(dof, tetrahedron), dp) * &
                    solution(global_dofs(dof, tetrahedron))
            end do
            curl_value = matmul(physical_curls, local_dofs)
            status = 0
            return
        end do
    end subroutine probe_box_curl

    pure subroutine box_reluctivity(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value(3, 3)

        value = 0.0_dp
        value(1, 1) = 1.0_dp
        value(2, 2) = 1.0_dp
        value(3, 3) = z
        associate (unused => [x, y])
            if (size(unused) /= 2) error stop
        end associate
    end subroutine box_reluctivity

    pure subroutine axial_source(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value(3)

        value = [0.0_dp, 0.0_dp, 1.0_dp]
        associate (unused => [x, y, z])
            if (size(unused) /= 3) error stop
        end associate
    end subroutine axial_source

end module fortfem_magnetic_box_3d
