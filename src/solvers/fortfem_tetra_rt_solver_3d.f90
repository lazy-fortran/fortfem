module fortfem_tetra_rt_solver_3d
    use fortfem_assembly_tetra_rt_arbitrary_order_3d, only: &
        assemble_tetra_rt_div_mass_csc, assemble_tetra_rt_vector_load
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_solve_csc
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
            mass_coefficient, solution, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree
        procedure(vector_source_3d) :: source
        real(dp), intent(in) :: divergence_coefficient, mass_coefficient
        real(dp), allocatable, intent(out) :: solution(:)
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: matrix
        real(dp), allocatable :: right_hand_side(:)
        integer :: solve_status

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Tetrahedral RT div-mass solve failed")
        if (degree < 0 .or. degree > 4) return
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
        call sparse_direct_solve_csc( &
            matrix%nrow, matrix%col_ptr, matrix%row_idx, matrix%val, &
            right_hand_side, solution, solve_status)
        if (solve_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INTERNAL_ERROR, &
                "Tetrahedral RT sparse solve failed")
            return
        end if
        call status_set(status, 0, "")
    end subroutine solve_tetra_rt_div_mass

end module fortfem_tetra_rt_solver_3d
