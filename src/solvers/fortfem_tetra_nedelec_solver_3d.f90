module fortfem_tetra_nedelec_solver_3d
    use fortfem_assembly_tetra_nedelec_3d, only: &
        assemble_tetra_nedelec_curl_mass_csc, &
        assemble_tetra_nedelec_vector_load_order
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_solve_csc
    use fortfem_tetra_nedelec_global_dof_map, only: &
        build_tetra_nedelec_dof_map
    use fortsparse, only: csc_t, FORTSPARSE_INTERNAL_ERROR, &
        FORTSPARSE_INVALID_MATRIX, fortsparse_status_t, status_set
    implicit none

    private

    public :: solve_tetra_nedelec_curl_mass

    abstract interface
        pure subroutine vector_source_3d(x, y, z, value)
            import :: dp
            real(dp), intent(in) :: x, y, z
            real(dp), intent(out) :: value(3)
        end subroutine vector_source_3d
    end interface

contains

    subroutine solve_tetra_nedelec_curl_mass( &
            vertices, tetrahedra, order, source, curl_coefficient, &
            mass_coefficient, solution, status, zero_tangential_boundary)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), order
        procedure(vector_source_3d) :: source
        real(dp), intent(in) :: curl_coefficient, mass_coefficient
        real(dp), allocatable, intent(out) :: solution(:)
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(in), optional :: zero_tangential_boundary

        type(csc_t) :: matrix
        logical, allocatable :: constrained(:)
        real(dp), allocatable :: right_hand_side(:)
        integer :: solve_status
        logical :: impose_boundary

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Tetrahedral Nedelec curl-mass solve failed")
        if (order < 1 .or. order > 4) return
        if (mass_coefficient <= 0.0_dp) return
        call assemble_tetra_nedelec_curl_mass_csc( &
            vertices, tetrahedra, matrix, status, curl_coefficient, &
            mass_coefficient, order)
        if (status%code /= 0) return
        call assemble_tetra_nedelec_vector_load_order( &
            vertices, tetrahedra, order, source, right_hand_side, status)
        if (status%code /= 0) return
        if (size(right_hand_side) /= matrix%nrow) return
        allocate(solution(matrix%nrow))
        impose_boundary = .false.
        if (present(zero_tangential_boundary)) then
            impose_boundary = zero_tangential_boundary
        end if
        if (impose_boundary) then
            call build_zero_tangential_mask( &
                tetrahedra, order, constrained, solve_status)
            if (solve_status /= 0) return
            call solve_zero_constrained( &
                matrix, right_hand_side, constrained, solution, solve_status)
        else
            call sparse_direct_solve_csc( &
                matrix%nrow, matrix%col_ptr, matrix%row_idx, matrix%val, &
                right_hand_side, solution, solve_status)
        end if
        if (solve_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INTERNAL_ERROR, &
                "Tetrahedral Nedelec sparse solve failed")
            return
        end if
        call status_set(status, 0, "")
    end subroutine solve_tetra_nedelec_curl_mass

    subroutine build_zero_tangential_mask( &
            tetrahedra, order, constrained, status)
        integer, intent(in) :: tetrahedra(:, :), order
        logical, allocatable, intent(out) :: constrained(:)
        integer, intent(out) :: status

        integer, allocatable :: edge_orientations(:, :), edges(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :)
        integer :: dof, edge, edge_dof_count, face, face_dof_count
        integer :: incidence, tetrahedron

        call build_tetra_nedelec_dof_map( &
            order, tetrahedra, edges, faces, global_dofs, edge_orientations, &
            face_permutations, status)
        if (status /= 0) return
        allocate(constrained(maxval(global_dofs)))
        constrained = .false.
        edge_dof_count = order*size(edges, 2)
        face_dof_count = order*(order - 1)
        do face = 1, size(faces, 2)
            incidence = 0
            do tetrahedron = 1, size(tetrahedra, 2)
                if (all_vertices_in_cell( &
                    faces(:, face), tetrahedra(:, tetrahedron))) then
                    incidence = incidence + 1
                end if
            end do
            if (incidence /= 1) cycle
            do edge = 1, size(edges, 2)
                if (.not. all_vertices_in_face( &
                    edges(:, edge), faces(:, face))) cycle
                do dof = 1, order
                    constrained((edge - 1)*order + dof) = .true.
                end do
            end do
            do dof = 1, face_dof_count
                constrained(edge_dof_count + &
                    (face - 1)*face_dof_count + dof) = .true.
            end do
        end do
        status = 0
    end subroutine build_zero_tangential_mask

    pure logical function all_vertices_in_cell(face, cell) result(found)
        integer, intent(in) :: face(3), cell(4)

        integer :: vertex

        found = .true.
        do vertex = 1, 3
            found = found .and. any(cell == face(vertex))
        end do
    end function all_vertices_in_cell

    pure logical function all_vertices_in_face(edge, face) result(found)
        integer, intent(in) :: edge(2), face(3)

        found = any(face == edge(1)) .and. any(face == edge(2))
    end function all_vertices_in_face

    subroutine solve_zero_constrained( &
            matrix, rhs, constrained, solution, status)
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: rhs(:)
        logical, intent(in) :: constrained(:)
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        integer, allocatable :: column_pointers(:), free_dofs(:)
        integer, allocatable :: free_index(:), row_indices(:)
        real(dp), allocatable :: reduced_rhs(:), reduced_solution(:)
        real(dp), allocatable :: values(:)
        integer :: column, entry, free_count, reduced_entry, row

        status = 1
        solution = 0.0_dp
        free_count = count(.not. constrained)
        if (free_count < 1) return
        allocate(free_dofs(free_count), free_index(matrix%nrow))
        free_count = 0
        do column = 1, matrix%ncol
            if (constrained(column)) cycle
            free_count = free_count + 1
            free_dofs(free_count) = column
        end do
        free_index = 0
        do column = 1, free_count
            free_index(free_dofs(column)) = column
        end do
        allocate(column_pointers(free_count + 1))
        reduced_entry = 0
        do column = 1, free_count
            column_pointers(column) = reduced_entry + 1
            do entry = matrix%col_ptr(free_dofs(column)), &
                    matrix%col_ptr(free_dofs(column) + 1) - 1
                if (.not. constrained(matrix%row_idx(entry))) then
                    reduced_entry = reduced_entry + 1
                end if
            end do
        end do
        column_pointers(free_count + 1) = reduced_entry + 1
        allocate(row_indices(reduced_entry), values(reduced_entry))
        reduced_entry = 0
        do column = 1, free_count
            do entry = matrix%col_ptr(free_dofs(column)), &
                    matrix%col_ptr(free_dofs(column) + 1) - 1
                row = matrix%row_idx(entry)
                if (constrained(row)) cycle
                reduced_entry = reduced_entry + 1
                row_indices(reduced_entry) = free_index(row)
                values(reduced_entry) = matrix%val(entry)
            end do
        end do
        allocate(reduced_rhs(free_count), reduced_solution(free_count))
        reduced_rhs = rhs(free_dofs)
        call sparse_direct_solve_csc( &
            free_count, column_pointers, row_indices, values, reduced_rhs, &
            reduced_solution, status)
        if (status /= 0) return
        solution(free_dofs) = reduced_solution
    end subroutine solve_zero_constrained

end module fortfem_tetra_nedelec_solver_3d
