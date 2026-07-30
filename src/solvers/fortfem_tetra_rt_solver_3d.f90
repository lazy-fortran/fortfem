module fortfem_tetra_rt_solver_3d
    use fortfem_assembly_tetra_rt_arbitrary_order_3d, only: &
        assemble_tetra_rt_div_mass_csc, assemble_tetra_rt_vector_load
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_solve_csc, &
        sparse_direct_solve_zero_constrained
    use fortfem_tetra_rt_global_dof_map, only: build_tetra_rt_dof_map
    use fortsparse, only: csc_t, FORTSPARSE_INTERNAL_ERROR, &
        FORTSPARSE_INVALID_MATRIX, fortsparse_status_t, status_set
    implicit none

    private

    public :: solve_tetra_rt_div_mass

    abstract interface
        pure subroutine vector_source_3d(x, y, z, value)
            import :: dp
            real(dp), intent(in) :: x, y, z
            real(dp), intent(out) :: value(3)
        end subroutine vector_source_3d
    end interface

contains

    subroutine solve_tetra_rt_div_mass( &
            vertices, tetrahedra, degree, source, divergence_coefficient, &
            mass_coefficient, solution, status, zero_normal_boundary)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree
        procedure(vector_source_3d) :: source
        real(dp), intent(in) :: divergence_coefficient, mass_coefficient
        real(dp), allocatable, intent(out) :: solution(:)
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(in), optional :: zero_normal_boundary

        type(csc_t) :: matrix
        logical, allocatable :: constrained(:)
        real(dp), allocatable :: right_hand_side(:)
        integer :: solve_status
        logical :: impose_boundary

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Tetrahedral RT div-mass solve failed")
        if (degree < 0) return
        if (mass_coefficient <= 0.0_dp) return
        call assemble_tetra_rt_div_mass_csc( &
            vertices, tetrahedra, degree, 2*degree + 4, matrix, status, &
            divergence_coefficient, mass_coefficient)
        if (status%code /= 0) return
        call assemble_tetra_rt_vector_load( &
            vertices, tetrahedra, degree, 2*degree + 4, source, &
            right_hand_side, status)
        if (status%code /= 0) return
        if (size(right_hand_side) /= matrix%nrow) return
        allocate(solution(matrix%nrow))
        impose_boundary = .false.
        if (present(zero_normal_boundary)) then
            impose_boundary = zero_normal_boundary
        end if
        if (impose_boundary) then
            call build_zero_normal_mask( &
                tetrahedra, degree, constrained, solve_status)
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
                "Tetrahedral RT sparse solve failed")
            return
        end if
        call status_set(status, 0, "")
    end subroutine solve_tetra_rt_div_mass

    subroutine build_zero_normal_mask( &
            tetrahedra, degree, constrained, status)
        integer, intent(in) :: tetrahedra(:, :), degree
        logical, allocatable, intent(out) :: constrained(:)
        integer, intent(out) :: status

        integer, allocatable :: face_orientations(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :)
        integer :: dof, face, face_dof_count, incidence, tetrahedron

        call build_tetra_rt_dof_map( &
            degree, tetrahedra, faces, global_dofs, face_orientations, &
            face_permutations, status)
        if (status /= 0) return
        allocate(constrained(maxval(global_dofs)))
        constrained = .false.
        face_dof_count = (degree + 1)*(degree + 2)/2
        do face = 1, size(faces, 2)
            incidence = 0
            do tetrahedron = 1, size(tetrahedra, 2)
                if (all_vertices_in_cell( &
                    faces(:, face), tetrahedra(:, tetrahedron))) then
                    incidence = incidence + 1
                end if
            end do
            if (incidence /= 1) cycle
            do dof = 1, face_dof_count
                constrained((face - 1)*face_dof_count + dof) = .true.
            end do
        end do
        status = 0
    end subroutine build_zero_normal_mask

    pure logical function all_vertices_in_cell(face, cell) result(found)
        integer, intent(in) :: face(3), cell(4)

        integer :: vertex

        found = .true.
        do vertex = 1, 3
            found = found .and. any(cell == face(vertex))
        end do
    end function all_vertices_in_cell

end module fortfem_tetra_rt_solver_3d
