module fortfem_tetra_nedelec_pml_state_3d
    !! Differentiable complex Nedelec PML state with optional boundary form.
    use fortfem_assembly_tetra_nedelec_3d, only: &
        assemble_tetra_nedelec_pml_csc, &
        assemble_tetra_nedelec_pml_csc_jvp, &
        assemble_tetra_nedelec_pml_csc_vjp, &
        assemble_tetra_nedelec_curvilinear_pml_csc, &
        assemble_tetra_nedelec_curvilinear_pml_csc_jvp, &
        assemble_tetra_nedelec_curvilinear_pml_csc_vjp
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: &
        sparse_direct_solve_constrained_jvp, &
        sparse_direct_solve_constrained_vjp
    use fortsparse, only: csc_from_triplet, csc_z_t, fortsparse_status_t, &
        status_set
    implicit none

    private

    public :: solve_tetra_nedelec_pml_jvp
    public :: solve_tetra_nedelec_pml_vjp
    public :: solve_tetra_nedelec_curvilinear_pml_jvp
    public :: solve_tetra_nedelec_curvilinear_pml_vjp

contains

    subroutine solve_tetra_nedelec_pml_jvp( &
            vertices, tetrahedra, order, stretch, wave_number, volume_load, &
            dirichlet_dofs, dirichlet_values, vertices_dot, stretch_dot, &
            wave_number_dot, volume_load_dot, dirichlet_values_dot, &
            solution_dot, status, boundary_operator_dofs, boundary_operator, &
            boundary_operator_dot)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        integer, intent(in) :: tetrahedra(:, :), order
        complex(dp), intent(in) :: stretch(:, :), stretch_dot(:, :)
        real(dp), intent(in) :: wave_number, wave_number_dot
        complex(dp), intent(in) :: volume_load(:), volume_load_dot(:)
        integer, intent(in) :: dirichlet_dofs(:)
        complex(dp), intent(in) :: dirichlet_values(:), dirichlet_values_dot(:)
        complex(dp), allocatable, intent(out) :: solution_dot(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: boundary_operator_dofs(:)
        complex(dp), intent(in), optional :: boundary_operator(:, :)
        complex(dp), intent(in), optional :: boundary_operator_dot(:, :)

        type(csc_z_t) :: pml_matrix, pml_matrix_dot
        type(csc_z_t) :: matrix, matrix_dot
        type(fortsparse_status_t) :: sparse_status
        logical, allocatable :: constrained(:)
        complex(dp), allocatable :: constrained_values(:)
        complex(dp), allocatable :: constrained_values_dot(:)
        logical :: has_boundary

        status = 1
        if (allocated(solution_dot)) deallocate(solution_dot)
        has_boundary = present(boundary_operator_dofs) .or. &
            present(boundary_operator) .or. present(boundary_operator_dot)
        if (has_boundary .and. (.not. present(boundary_operator_dofs) .or. &
            .not. present(boundary_operator) .or. &
            .not. present(boundary_operator_dot))) return
        call assemble_tetra_nedelec_pml_csc( &
            vertices, tetrahedra, order, stretch, wave_number, pml_matrix, &
            sparse_status)
        if (sparse_status%code /= 0) return
        call assemble_tetra_nedelec_pml_csc_jvp( &
            vertices, tetrahedra, order, stretch, wave_number, vertices_dot, &
            stretch_dot, wave_number_dot, pml_matrix_dot, sparse_status)
        if (sparse_status%code /= 0) return
        call combine_boundary_form( &
            pml_matrix, boundary_operator_dofs, boundary_operator, matrix, &
            sparse_status)
        if (sparse_status%code /= 0) return
        call combine_boundary_form( &
            pml_matrix_dot, boundary_operator_dofs, boundary_operator_dot, &
            matrix_dot, sparse_status)
        if (sparse_status%code /= 0) return
        if (.not. same_csc_pattern(matrix, matrix_dot)) return
        if (size(volume_load) /= matrix%nrow .or. &
            size(volume_load_dot) /= matrix%nrow) return
        call build_constraints( &
            matrix%nrow, dirichlet_dofs, dirichlet_values, &
            dirichlet_values_dot, constrained, constrained_values, &
            constrained_values_dot, status)
        if (status /= 0) return
        allocate(solution_dot(matrix%nrow))
        call sparse_direct_solve_constrained_jvp( &
            matrix, volume_load, constrained, constrained_values, &
            matrix_dot%val, volume_load_dot, constrained_values_dot, &
            solution_dot, status)
    end subroutine solve_tetra_nedelec_pml_jvp

    subroutine solve_tetra_nedelec_pml_vjp( &
            vertices, tetrahedra, order, stretch, wave_number, volume_load, &
            dirichlet_dofs, dirichlet_values, solution, solution_bar, &
            vertices_bar, stretch_bar, wave_number_bar, volume_load_bar, &
            dirichlet_values_bar, status, boundary_operator_dofs, &
            boundary_operator, boundary_operator_bar)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), order
        complex(dp), intent(in) :: stretch(:, :), volume_load(:)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: dirichlet_dofs(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), intent(in) :: solution(:), solution_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :), wave_number_bar
        complex(dp), intent(out) :: stretch_bar(:, :), volume_load_bar(:)
        complex(dp), intent(out) :: dirichlet_values_bar(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: boundary_operator_dofs(:)
        complex(dp), intent(in), optional :: boundary_operator(:, :)
        complex(dp), intent(out), optional :: boundary_operator_bar(:, :)

        type(csc_z_t) :: pml_matrix, matrix
        type(fortsparse_status_t) :: sparse_status
        logical, allocatable :: constrained(:)
        complex(dp), allocatable :: constrained_values(:)
        complex(dp), allocatable :: constrained_values_bar(:)
        complex(dp), allocatable :: matrix_values_bar(:), pml_values_bar(:)
        complex(dp), allocatable :: rhs_bar(:)
        logical :: has_boundary
        integer :: constraint

        vertices_bar = 0.0_dp
        stretch_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wave_number_bar = 0.0_dp
        volume_load_bar = cmplx(0.0_dp, 0.0_dp, dp)
        dirichlet_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(boundary_operator_bar)) then
            boundary_operator_bar = cmplx(0.0_dp, 0.0_dp, dp)
        end if
        status = 1
        has_boundary = present(boundary_operator_dofs) .or. &
            present(boundary_operator) .or. present(boundary_operator_bar)
        if (has_boundary .and. (.not. present(boundary_operator_dofs) .or. &
            .not. present(boundary_operator) .or. &
            .not. present(boundary_operator_bar))) return
        call assemble_tetra_nedelec_pml_csc( &
            vertices, tetrahedra, order, stretch, wave_number, pml_matrix, &
            sparse_status)
        if (sparse_status%code /= 0) return
        call combine_boundary_form( &
            pml_matrix, boundary_operator_dofs, boundary_operator, matrix, &
            sparse_status)
        if (sparse_status%code /= 0) return
        if (size(volume_load) /= matrix%nrow .or. &
            size(solution) /= matrix%nrow .or. &
            size(solution_bar) /= matrix%nrow .or. &
            size(volume_load_bar) /= matrix%nrow .or. &
            size(dirichlet_values_bar) /= size(dirichlet_values)) return
        call build_constraints( &
            matrix%nrow, dirichlet_dofs, dirichlet_values, &
            dirichlet_values, constrained, constrained_values, &
            constrained_values_bar, status)
        if (status /= 0) return
        allocate(matrix_values_bar(matrix%nnz), rhs_bar(matrix%nrow))
        call sparse_direct_solve_constrained_vjp( &
            matrix, volume_load, constrained, constrained_values, solution, &
            solution_bar, matrix_values_bar, rhs_bar, &
            constrained_values_bar, status)
        if (status /= 0) return
        volume_load_bar = rhs_bar
        do constraint = 1, size(dirichlet_dofs)
            dirichlet_values_bar(constraint) = &
                constrained_values_bar(dirichlet_dofs(constraint))
        end do
        allocate(pml_values_bar(pml_matrix%nnz))
        call restrict_matrix_bar( &
            pml_matrix, matrix, matrix_values_bar, pml_values_bar)
        call assemble_tetra_nedelec_pml_csc_vjp( &
            vertices, tetrahedra, order, stretch, wave_number, &
            pml_values_bar, vertices_bar, stretch_bar, wave_number_bar, &
            sparse_status)
        if (sparse_status%code /= 0) then
            status = sparse_status%code
            return
        end if
        if (has_boundary) call boundary_form_vjp( &
            matrix, matrix_values_bar, boundary_operator_dofs, &
            boundary_operator_bar)
        status = 0
    end subroutine solve_tetra_nedelec_pml_vjp

    subroutine solve_tetra_nedelec_curvilinear_pml_jvp( &
            vertices, tetrahedra, order, stretch, wave_number, volume_load, &
            dirichlet_dofs, dirichlet_values, vertices_dot, stretch_dot, &
            wave_number_dot, volume_load_dot, dirichlet_values_dot, &
            solution_dot, status, boundary_operator_dofs, boundary_operator, &
            boundary_operator_dot)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        integer, intent(in) :: tetrahedra(:, :), order
        complex(dp), intent(in) :: stretch(:, :, :), stretch_dot(:, :, :)
        real(dp), intent(in) :: wave_number, wave_number_dot
        complex(dp), intent(in) :: volume_load(:), volume_load_dot(:)
        integer, intent(in) :: dirichlet_dofs(:)
        complex(dp), intent(in) :: dirichlet_values(:), dirichlet_values_dot(:)
        complex(dp), allocatable, intent(out) :: solution_dot(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: boundary_operator_dofs(:)
        complex(dp), intent(in), optional :: boundary_operator(:, :)
        complex(dp), intent(in), optional :: boundary_operator_dot(:, :)

        type(csc_z_t) :: pml_matrix, pml_matrix_dot
        type(csc_z_t) :: matrix, matrix_dot
        type(fortsparse_status_t) :: sparse_status
        logical, allocatable :: constrained(:)
        complex(dp), allocatable :: constrained_values(:)
        complex(dp), allocatable :: constrained_values_dot(:)
        logical :: has_boundary

        status = 1
        if (allocated(solution_dot)) deallocate(solution_dot)
        has_boundary = present(boundary_operator_dofs) .or. &
            present(boundary_operator) .or. present(boundary_operator_dot)
        if (has_boundary .and. (.not. present(boundary_operator_dofs) .or. &
            .not. present(boundary_operator) .or. &
            .not. present(boundary_operator_dot))) return
        call assemble_tetra_nedelec_curvilinear_pml_csc( &
            vertices, tetrahedra, order, stretch, wave_number, pml_matrix, &
            sparse_status)
        if (sparse_status%code /= 0) return
        call assemble_tetra_nedelec_curvilinear_pml_csc_jvp( &
            vertices, tetrahedra, order, stretch, wave_number, vertices_dot, &
            stretch_dot, wave_number_dot, pml_matrix_dot, sparse_status)
        if (sparse_status%code /= 0) return
        call combine_boundary_form( &
            pml_matrix, boundary_operator_dofs, boundary_operator, matrix, &
            sparse_status)
        if (sparse_status%code /= 0) return
        call combine_boundary_form( &
            pml_matrix_dot, boundary_operator_dofs, boundary_operator_dot, &
            matrix_dot, sparse_status)
        if (sparse_status%code /= 0) return
        if (.not. same_csc_pattern(matrix, matrix_dot)) return
        if (size(volume_load) /= matrix%nrow .or. &
            size(volume_load_dot) /= matrix%nrow) return
        call build_constraints( &
            matrix%nrow, dirichlet_dofs, dirichlet_values, &
            dirichlet_values_dot, constrained, constrained_values, &
            constrained_values_dot, status)
        if (status /= 0) return
        allocate(solution_dot(matrix%nrow))
        call sparse_direct_solve_constrained_jvp( &
            matrix, volume_load, constrained, constrained_values, &
            matrix_dot%val, volume_load_dot, constrained_values_dot, &
            solution_dot, status)
    end subroutine solve_tetra_nedelec_curvilinear_pml_jvp

    subroutine solve_tetra_nedelec_curvilinear_pml_vjp( &
            vertices, tetrahedra, order, stretch, wave_number, volume_load, &
            dirichlet_dofs, dirichlet_values, solution, solution_bar, &
            vertices_bar, stretch_bar, wave_number_bar, volume_load_bar, &
            dirichlet_values_bar, status, boundary_operator_dofs, &
            boundary_operator, boundary_operator_bar)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), order
        complex(dp), intent(in) :: stretch(:, :, :), volume_load(:)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: dirichlet_dofs(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), intent(in) :: solution(:), solution_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :), wave_number_bar
        complex(dp), intent(out) :: stretch_bar(:, :, :), volume_load_bar(:)
        complex(dp), intent(out) :: dirichlet_values_bar(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: boundary_operator_dofs(:)
        complex(dp), intent(in), optional :: boundary_operator(:, :)
        complex(dp), intent(out), optional :: boundary_operator_bar(:, :)

        type(csc_z_t) :: pml_matrix, matrix
        type(fortsparse_status_t) :: sparse_status
        logical, allocatable :: constrained(:)
        complex(dp), allocatable :: constrained_values(:)
        complex(dp), allocatable :: constrained_values_bar(:)
        complex(dp), allocatable :: matrix_values_bar(:), pml_values_bar(:)
        complex(dp), allocatable :: rhs_bar(:)
        logical :: has_boundary
        integer :: constraint

        vertices_bar = 0.0_dp
        stretch_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wave_number_bar = 0.0_dp
        volume_load_bar = cmplx(0.0_dp, 0.0_dp, dp)
        dirichlet_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        if (present(boundary_operator_bar)) then
            boundary_operator_bar = cmplx(0.0_dp, 0.0_dp, dp)
        end if
        status = 1
        has_boundary = present(boundary_operator_dofs) .or. &
            present(boundary_operator) .or. present(boundary_operator_bar)
        if (has_boundary .and. (.not. present(boundary_operator_dofs) .or. &
            .not. present(boundary_operator) .or. &
            .not. present(boundary_operator_bar))) return
        call assemble_tetra_nedelec_curvilinear_pml_csc( &
            vertices, tetrahedra, order, stretch, wave_number, pml_matrix, &
            sparse_status)
        if (sparse_status%code /= 0) return
        call combine_boundary_form( &
            pml_matrix, boundary_operator_dofs, boundary_operator, matrix, &
            sparse_status)
        if (sparse_status%code /= 0) return
        if (size(volume_load) /= matrix%nrow .or. &
            size(solution) /= matrix%nrow .or. &
            size(solution_bar) /= matrix%nrow .or. &
            size(volume_load_bar) /= matrix%nrow .or. &
            size(dirichlet_values_bar) /= size(dirichlet_values)) return
        call build_constraints( &
            matrix%nrow, dirichlet_dofs, dirichlet_values, &
            dirichlet_values, constrained, constrained_values, &
            constrained_values_bar, status)
        if (status /= 0) return
        allocate(matrix_values_bar(matrix%nnz), rhs_bar(matrix%nrow))
        call sparse_direct_solve_constrained_vjp( &
            matrix, volume_load, constrained, constrained_values, solution, &
            solution_bar, matrix_values_bar, rhs_bar, &
            constrained_values_bar, status)
        if (status /= 0) return
        volume_load_bar = rhs_bar
        do constraint = 1, size(dirichlet_dofs)
            dirichlet_values_bar(constraint) = &
                constrained_values_bar(dirichlet_dofs(constraint))
        end do
        allocate(pml_values_bar(pml_matrix%nnz))
        call restrict_matrix_bar( &
            pml_matrix, matrix, matrix_values_bar, pml_values_bar)
        call assemble_tetra_nedelec_curvilinear_pml_csc_vjp( &
            vertices, tetrahedra, order, stretch, wave_number, &
            pml_values_bar, vertices_bar, stretch_bar, wave_number_bar, &
            sparse_status)
        if (sparse_status%code /= 0) then
            status = sparse_status%code
            return
        end if
        if (has_boundary) call boundary_form_vjp( &
            matrix, matrix_values_bar, boundary_operator_dofs, &
            boundary_operator_bar)
        status = 0
    end subroutine solve_tetra_nedelec_curvilinear_pml_vjp

    subroutine combine_boundary_form( &
            pml_matrix, boundary_operator_dofs, boundary_operator, matrix, &
            status)
        type(csc_z_t), intent(in) :: pml_matrix
        integer, intent(in), optional :: boundary_operator_dofs(:)
        complex(dp), intent(in), optional :: boundary_operator(:, :)
        type(csc_z_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: rows(:), columns(:)
        complex(dp), allocatable :: values(:)
        integer :: column, entry, operator_entry, row

        if (.not. present(boundary_operator_dofs) .and. &
            .not. present(boundary_operator)) then
            matrix = pml_matrix
            call status_set(status, 0, "")
            return
        end if
        if (.not. present(boundary_operator_dofs) .or. &
            .not. present(boundary_operator)) then
            call status_set(status, 1, "Boundary form arguments differ")
            return
        end if
        if (size(boundary_operator, 1) /= size(boundary_operator_dofs) .or. &
            size(boundary_operator, 2) /= size(boundary_operator_dofs) .or. &
            any(boundary_operator_dofs < 1) .or. &
            any(boundary_operator_dofs > pml_matrix%nrow)) then
            call status_set(status, 1, "Boundary form shape is invalid")
            return
        end if
        allocate(rows(pml_matrix%nnz + size(boundary_operator_dofs)**2))
        allocate(columns(size(rows)), values(size(rows)))
        operator_entry = 0
        do column = 1, pml_matrix%ncol
            do entry = pml_matrix%col_ptr(column), &
                    pml_matrix%col_ptr(column + 1) - 1
                operator_entry = operator_entry + 1
                rows(operator_entry) = pml_matrix%row_idx(entry)
                columns(operator_entry) = column
                values(operator_entry) = pml_matrix%val(entry)
            end do
        end do
        do column = 1, size(boundary_operator_dofs)
            do row = 1, size(boundary_operator_dofs)
                operator_entry = operator_entry + 1
                rows(operator_entry) = boundary_operator_dofs(row)
                columns(operator_entry) = boundary_operator_dofs(column)
                values(operator_entry) = boundary_operator(row, column)
            end do
        end do
        call csc_from_triplet( &
            pml_matrix%nrow, pml_matrix%ncol, rows, columns, values, matrix, &
            status)
    end subroutine combine_boundary_form

    subroutine build_constraints( &
            n, dirichlet_dofs, dirichlet_values, dirichlet_values_dot, &
            constrained, values, values_dot, status)
        integer, intent(in) :: n, dirichlet_dofs(:)
        complex(dp), intent(in) :: dirichlet_values(:), dirichlet_values_dot(:)
        logical, allocatable, intent(out) :: constrained(:)
        complex(dp), allocatable, intent(out) :: values(:), values_dot(:)
        integer, intent(out) :: status

        integer :: constraint, dof

        status = 1
        if (size(dirichlet_values) /= size(dirichlet_dofs) .or. &
            size(dirichlet_values_dot) /= size(dirichlet_dofs)) return
        if (any(dirichlet_dofs < 1) .or. any(dirichlet_dofs > n)) return
        allocate(constrained(n), values(n), values_dot(n))
        constrained = .false.
        values = cmplx(0.0_dp, 0.0_dp, dp)
        values_dot = cmplx(0.0_dp, 0.0_dp, dp)
        do constraint = 1, size(dirichlet_dofs)
            dof = dirichlet_dofs(constraint)
            if (constrained(dof)) return
            constrained(dof) = .true.
            values(dof) = dirichlet_values(constraint)
            values_dot(dof) = dirichlet_values_dot(constraint)
        end do
        status = 0
    end subroutine build_constraints

    pure logical function same_csc_pattern(first, second) result(same)
        type(csc_z_t), intent(in) :: first, second

        same = first%nrow == second%nrow .and. first%ncol == second%ncol
        if (.not. same) return
        same = first%nnz == second%nnz
        if (.not. same) return
        same = all(first%col_ptr == second%col_ptr) .and. &
            all(first%row_idx == second%row_idx)
    end function same_csc_pattern

    subroutine restrict_matrix_bar( &
            pml_matrix, matrix, matrix_values_bar, pml_values_bar)
        type(csc_z_t), intent(in) :: pml_matrix, matrix
        complex(dp), intent(in) :: matrix_values_bar(:)
        complex(dp), intent(out) :: pml_values_bar(:)

        integer :: column, entry, global_entry, row

        pml_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        do column = 1, pml_matrix%ncol
            do entry = pml_matrix%col_ptr(column), &
                    pml_matrix%col_ptr(column + 1) - 1
                row = pml_matrix%row_idx(entry)
                global_entry = csc_entry(matrix, row, column)
                if (global_entry > 0) pml_values_bar(entry) = &
                    matrix_values_bar(global_entry)
            end do
        end do
    end subroutine restrict_matrix_bar

    subroutine boundary_form_vjp( &
            matrix, matrix_values_bar, boundary_operator_dofs, &
            boundary_operator_bar)
        type(csc_z_t), intent(in) :: matrix
        complex(dp), intent(in) :: matrix_values_bar(:)
        integer, intent(in) :: boundary_operator_dofs(:)
        complex(dp), intent(out) :: boundary_operator_bar(:, :)

        integer :: column, global_entry, row

        boundary_operator_bar = cmplx(0.0_dp, 0.0_dp, dp)
        do column = 1, size(boundary_operator_dofs)
            do row = 1, size(boundary_operator_dofs)
                global_entry = csc_entry( &
                    matrix, boundary_operator_dofs(row), &
                    boundary_operator_dofs(column))
                if (global_entry > 0) boundary_operator_bar(row, column) = &
                    matrix_values_bar(global_entry)
            end do
        end do
    end subroutine boundary_form_vjp

    pure integer function csc_entry(matrix, row, column) result(entry_index)
        type(csc_z_t), intent(in) :: matrix
        integer, intent(in) :: row, column

        integer :: entry

        entry_index = 0
        if (column < 1 .or. column > matrix%ncol) return
        do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
            if (matrix%row_idx(entry) == row) then
                entry_index = entry
                return
            end if
        end do
    end function csc_entry

end module fortfem_tetra_nedelec_pml_state_3d
