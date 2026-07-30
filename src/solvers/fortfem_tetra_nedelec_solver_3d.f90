module fortfem_tetra_nedelec_solver_3d
    use fortfem_assembly_tetra_nedelec_3d, only: &
        assemble_tetra_nedelec_curl_mass_csc, &
        assemble_tetra_nedelec_vector_load_order
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_solve_csc, &
        sparse_direct_solve_zero_constrained
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
            call sparse_direct_solve_zero_constrained( &
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

end module fortfem_tetra_nedelec_solver_3d
