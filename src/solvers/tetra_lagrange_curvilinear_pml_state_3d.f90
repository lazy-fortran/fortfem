module fortfem_tetra_lagrange_curvilinear_pml_state_3d
    !! Differentiable constrained scalar curvilinear PML state.
    use fortfem_assembly_tetra_lagrange_curvilinear_pml_3d, only: &
        assemble_tetra_lagrange_curvilinear_pml_csc, &
        assemble_tetra_lagrange_curvilinear_pml_csc_jvp, &
        assemble_tetra_lagrange_curvilinear_pml_csc_vjp
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: &
        sparse_direct_solve_constrained, &
        sparse_direct_solve_constrained_jvp, &
        sparse_direct_solve_constrained_vjp
    use fortsparse, only: csc_z_t, fortsparse_status_t
    implicit none
    private

    public :: solve_tetra_lagrange_curvilinear_pml
    public :: solve_tetra_lagrange_curvilinear_pml_jvp
    public :: solve_tetra_lagrange_curvilinear_pml_vjp

contains

    subroutine solve_tetra_lagrange_curvilinear_pml( &
            vertices, tetrahedra, degree, stretch, wave_number, load, &
            constrained, constrained_values, solution, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree
        complex(dp), intent(in) :: stretch(:, :, :), load(:)
        real(dp), intent(in) :: wave_number
        logical, intent(in) :: constrained(:)
        complex(dp), intent(in) :: constrained_values(:)
        complex(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: status

        type(csc_z_t) :: matrix
        type(fortsparse_status_t) :: sparse_status

        status = 1
        if (allocated(solution)) deallocate(solution)
        call assemble_tetra_lagrange_curvilinear_pml_csc( &
            vertices, tetrahedra, degree, stretch, wave_number, matrix, &
            sparse_status)
        if (sparse_status%code /= 0) return
        if (.not. valid_state_shapes( &
            matrix, load, constrained, constrained_values)) return
        allocate(solution(matrix%nrow))
        call sparse_direct_solve_constrained( &
            matrix, load, constrained, constrained_values, solution, status)
    end subroutine solve_tetra_lagrange_curvilinear_pml

    subroutine solve_tetra_lagrange_curvilinear_pml_jvp( &
            vertices, tetrahedra, degree, stretch, wave_number, load, &
            constrained, constrained_values, vertices_dot, stretch_dot, &
            wave_number_dot, load_dot, constrained_values_dot, solution_dot, &
            status)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree
        complex(dp), intent(in) :: stretch(:, :, :), stretch_dot(:, :, :)
        real(dp), intent(in) :: wave_number, wave_number_dot
        complex(dp), intent(in) :: load(:), load_dot(:)
        logical, intent(in) :: constrained(:)
        complex(dp), intent(in) :: constrained_values(:)
        complex(dp), intent(in) :: constrained_values_dot(:)
        complex(dp), allocatable, intent(out) :: solution_dot(:)
        integer, intent(out) :: status

        type(csc_z_t) :: matrix, matrix_dot
        type(fortsparse_status_t) :: sparse_status

        status = 1
        if (allocated(solution_dot)) deallocate(solution_dot)
        call assemble_tetra_lagrange_curvilinear_pml_csc( &
            vertices, tetrahedra, degree, stretch, wave_number, matrix, &
            sparse_status)
        if (sparse_status%code /= 0) return
        call assemble_tetra_lagrange_curvilinear_pml_csc_jvp( &
            vertices, tetrahedra, degree, stretch, wave_number, vertices_dot, &
            stretch_dot, wave_number_dot, matrix_dot, sparse_status)
        if (sparse_status%code /= 0) return
        if (.not. same_csc_pattern(matrix, matrix_dot)) return
        if (.not. valid_state_shapes( &
            matrix, load, constrained, constrained_values)) return
        if (size(vertices_dot, 1) /= size(vertices, 1) .or. &
            size(vertices_dot, 2) /= size(vertices, 2)) return
        if (any(shape(stretch_dot) /= shape(stretch))) return
        if (size(load_dot) /= size(load)) return
        if (size(constrained_values_dot) /= size(constrained_values)) return
        allocate(solution_dot(matrix%nrow))
        call sparse_direct_solve_constrained_jvp( &
            matrix, load, constrained, constrained_values, matrix_dot%val, &
            load_dot, constrained_values_dot, solution_dot, status)
    end subroutine solve_tetra_lagrange_curvilinear_pml_jvp

    subroutine solve_tetra_lagrange_curvilinear_pml_vjp( &
            vertices, tetrahedra, degree, stretch, wave_number, load, &
            constrained, constrained_values, solution, solution_bar, &
            vertices_bar, stretch_bar, wave_number_bar, load_bar, &
            constrained_values_bar, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree
        complex(dp), intent(in) :: stretch(:, :, :), load(:)
        real(dp), intent(in) :: wave_number
        logical, intent(in) :: constrained(:)
        complex(dp), intent(in) :: constrained_values(:), solution(:)
        complex(dp), intent(in) :: solution_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :), wave_number_bar
        complex(dp), intent(out) :: stretch_bar(:, :, :), load_bar(:)
        complex(dp), intent(out) :: constrained_values_bar(:)
        integer, intent(out) :: status

        type(csc_z_t) :: matrix
        type(fortsparse_status_t) :: sparse_status
        complex(dp), allocatable :: matrix_values_bar(:), rhs_bar(:)
        complex(dp), allocatable :: constrained_values_bar_local(:)

        vertices_bar = 0.0_dp
        stretch_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wave_number_bar = 0.0_dp
        load_bar = cmplx(0.0_dp, 0.0_dp, dp)
        constrained_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        call assemble_tetra_lagrange_curvilinear_pml_csc( &
            vertices, tetrahedra, degree, stretch, wave_number, matrix, &
            sparse_status)
        if (sparse_status%code /= 0) return
        if (.not. valid_state_shapes( &
            matrix, load, constrained, constrained_values)) return
        if (size(solution) /= matrix%nrow .or. &
            size(solution_bar) /= matrix%nrow) return
        if (any(shape(vertices_bar) /= shape(vertices))) return
        if (any(shape(stretch_bar) /= shape(stretch))) return
        if (size(load_bar) /= size(load)) return
        if (size(constrained_values_bar) /= size(constrained_values)) return
        allocate(matrix_values_bar(matrix%nnz), rhs_bar(matrix%nrow))
        allocate(constrained_values_bar_local(matrix%nrow))
        call sparse_direct_solve_constrained_vjp( &
            matrix, load, constrained, constrained_values, solution, &
            solution_bar, matrix_values_bar, rhs_bar, &
            constrained_values_bar_local, status)
        if (status /= 0) return
        call assemble_tetra_lagrange_curvilinear_pml_csc_vjp( &
            vertices, tetrahedra, degree, stretch, wave_number, &
            matrix_values_bar, vertices_bar, stretch_bar, wave_number_bar, &
            sparse_status)
        if (sparse_status%code /= 0) then
            status = sparse_status%code
            return
        end if
        load_bar = rhs_bar
        constrained_values_bar = constrained_values_bar_local
        status = 0
    end subroutine solve_tetra_lagrange_curvilinear_pml_vjp

    pure logical function valid_state_shapes( &
            matrix, load, constrained, constrained_values) result(valid)
        type(csc_z_t), intent(in) :: matrix
        complex(dp), intent(in) :: load(:), constrained_values(:)
        logical, intent(in) :: constrained(:)

        valid = matrix%nrow == matrix%ncol .and. matrix%nrow > 0 .and. &
            size(load) == matrix%nrow .and. &
            size(constrained) == matrix%nrow .and. &
            size(constrained_values) == matrix%nrow
    end function valid_state_shapes

    pure logical function same_csc_pattern(left, right) result(same)
        type(csc_z_t), intent(in) :: left, right

        same = left%nrow == right%nrow .and. left%ncol == right%ncol .and. &
            left%nnz == right%nnz
        if (.not. same) return
        same = all(left%col_ptr == right%col_ptr) .and. &
            all(left%row_idx == right%row_idx)
    end function same_csc_pattern

end module fortfem_tetra_lagrange_curvilinear_pml_state_3d
