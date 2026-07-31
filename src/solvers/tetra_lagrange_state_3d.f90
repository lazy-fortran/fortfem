module fortfem_tetra_lagrange_state_3d
    !! Differentiable constrained tetrahedral H1 diffusion-reaction state.
    use fortfem_assembly_tetra_lagrange_arbitrary_order_3d, only: &
        assemble_tetra_lagrange_stiffness_csc, &
        assemble_tetra_lagrange_stiffness_csc_jvp, &
        assemble_tetra_lagrange_stiffness_csc_vjp
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_solve_constrained, &
        sparse_direct_solve_constrained_jvp, &
        sparse_direct_solve_constrained_vjp
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none

    private

    public :: solve_tetra_lagrange_state
    public :: solve_tetra_lagrange_state_jvp
    public :: solve_tetra_lagrange_state_vjp

contains

    subroutine solve_tetra_lagrange_state( &
            vertices, tetrahedra, degree, quadrature_degree, &
            stiffness_coefficient, mass_coefficient, load, constrained, &
            constrained_values, solution, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        real(dp), intent(in) :: stiffness_coefficient, mass_coefficient
        real(dp), intent(in) :: load(:), constrained_values(:)
        logical, intent(in) :: constrained(:)
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        type(csc_t) :: matrix
        type(fortsparse_status_t) :: sparse_status

        solution = 0.0_dp
        status = 1
        if (.not. valid_coefficients( &
            stiffness_coefficient, mass_coefficient)) return
        call assemble_tetra_lagrange_stiffness_csc( &
            vertices, tetrahedra, degree, quadrature_degree, matrix, &
            sparse_status, stiffness_coefficient, mass_coefficient)
        if (sparse_status%code /= 0) return
        if (.not. valid_state_shapes( &
            matrix, load, constrained, constrained_values, solution)) return
        call sparse_direct_solve_constrained( &
            matrix, load, constrained, constrained_values, solution, status)
    end subroutine solve_tetra_lagrange_state

    subroutine solve_tetra_lagrange_state_jvp( &
            vertices, tetrahedra, degree, quadrature_degree, &
            stiffness_coefficient, mass_coefficient, load, constrained, &
            constrained_values, vertices_dot, stiffness_coefficient_dot, &
            mass_coefficient_dot, load_dot, constrained_values_dot, &
            solution_dot, status)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        real(dp), intent(in) :: stiffness_coefficient, mass_coefficient
        real(dp), intent(in) :: stiffness_coefficient_dot
        real(dp), intent(in) :: mass_coefficient_dot
        real(dp), intent(in) :: load(:), load_dot(:)
        real(dp), intent(in) :: constrained_values(:)
        real(dp), intent(in) :: constrained_values_dot(:)
        logical, intent(in) :: constrained(:)
        real(dp), intent(out) :: solution_dot(:)
        integer, intent(out) :: status

        type(csc_t) :: matrix, matrix_dot
        type(fortsparse_status_t) :: sparse_status

        solution_dot = 0.0_dp
        status = 1
        if (.not. valid_coefficients( &
            stiffness_coefficient, mass_coefficient)) return
        call assemble_tetra_lagrange_stiffness_csc( &
            vertices, tetrahedra, degree, quadrature_degree, matrix, &
            sparse_status, stiffness_coefficient, mass_coefficient)
        if (sparse_status%code /= 0) return
        call assemble_tetra_lagrange_stiffness_csc_jvp( &
            vertices, tetrahedra, degree, quadrature_degree, &
            stiffness_coefficient, mass_coefficient, vertices_dot, &
            stiffness_coefficient_dot, mass_coefficient_dot, matrix_dot, &
            sparse_status)
        if (sparse_status%code /= 0) return
        if (.not. valid_state_shapes( &
            matrix, load, constrained, constrained_values, solution_dot)) &
            return
        if (size(load_dot) /= size(load)) return
        if (size(constrained_values_dot) /= size(constrained_values)) return
        call sparse_direct_solve_constrained_jvp( &
            matrix, load, constrained, constrained_values, matrix_dot%val, &
            load_dot, constrained_values_dot, solution_dot, status)
    end subroutine solve_tetra_lagrange_state_jvp

    subroutine solve_tetra_lagrange_state_vjp( &
            vertices, tetrahedra, degree, quadrature_degree, &
            stiffness_coefficient, mass_coefficient, load, constrained, &
            constrained_values, solution, solution_bar, vertices_bar, &
            stiffness_coefficient_bar, mass_coefficient_bar, load_bar, &
            constrained_values_bar, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        real(dp), intent(in) :: stiffness_coefficient, mass_coefficient
        real(dp), intent(in) :: load(:), constrained_values(:)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: solution(:), solution_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :)
        real(dp), intent(out) :: stiffness_coefficient_bar
        real(dp), intent(out) :: mass_coefficient_bar
        real(dp), intent(out) :: load_bar(:), constrained_values_bar(:)
        integer, intent(out) :: status

        type(csc_t) :: matrix
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: matrix_values_bar(:)

        vertices_bar = 0.0_dp
        stiffness_coefficient_bar = 0.0_dp
        mass_coefficient_bar = 0.0_dp
        load_bar = 0.0_dp
        constrained_values_bar = 0.0_dp
        status = 1
        if (.not. valid_coefficients( &
            stiffness_coefficient, mass_coefficient)) return
        if (any(shape(vertices_bar) /= shape(vertices))) return
        call assemble_tetra_lagrange_stiffness_csc( &
            vertices, tetrahedra, degree, quadrature_degree, matrix, &
            sparse_status, stiffness_coefficient, mass_coefficient)
        if (sparse_status%code /= 0) return
        if (.not. valid_state_shapes( &
            matrix, load, constrained, constrained_values, solution)) return
        if (size(solution_bar) /= size(solution)) return
        if (size(load_bar) /= size(load)) return
        if (size(constrained_values_bar) /= size(constrained_values)) return
        allocate(matrix_values_bar(matrix%nnz))
        call sparse_direct_solve_constrained_vjp( &
            matrix, load, constrained, constrained_values, solution, &
            solution_bar, matrix_values_bar, load_bar, &
            constrained_values_bar, status)
        if (status /= 0) return
        call assemble_tetra_lagrange_stiffness_csc_vjp( &
            vertices, tetrahedra, degree, quadrature_degree, &
            stiffness_coefficient, mass_coefficient, matrix_values_bar, &
            vertices_bar, stiffness_coefficient_bar, mass_coefficient_bar, &
            sparse_status)
        status = sparse_status%code
    end subroutine solve_tetra_lagrange_state_vjp

    pure logical function valid_coefficients( &
            stiffness_coefficient, mass_coefficient) result(valid)
        real(dp), intent(in) :: stiffness_coefficient, mass_coefficient

        valid = stiffness_coefficient >= 0.0_dp .and. &
            mass_coefficient >= 0.0_dp
        if (.not. valid) return
        valid = stiffness_coefficient > 0.0_dp .or. mass_coefficient > 0.0_dp
    end function valid_coefficients

    pure logical function valid_state_shapes( &
            matrix, load, constrained, constrained_values, solution) &
            result(valid)
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: load(:), constrained_values(:), solution(:)
        logical, intent(in) :: constrained(:)

        valid = size(load) == matrix%nrow
        if (.not. valid) return
        valid = size(constrained) == matrix%nrow
        if (.not. valid) return
        valid = size(constrained_values) == matrix%nrow
        if (.not. valid) return
        valid = size(solution) == matrix%nrow
    end function valid_state_shapes

end module fortfem_tetra_lagrange_state_3d
