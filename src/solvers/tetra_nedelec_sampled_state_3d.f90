module fortfem_tetra_nedelec_sampled_state_3d
    !! Differentiable constrained curl--curl state with sampled forcing.
    use fortfem_assembly_tetra_nedelec_3d, only: &
        assemble_tetra_nedelec_curl_mass_csc, &
        assemble_tetra_nedelec_curl_mass_csc_jvp, &
        assemble_tetra_nedelec_curl_mass_csc_vjp, &
        assemble_tetra_nedelec_vector_load_samples, &
        assemble_tetra_nedelec_vector_load_samples_jvp, &
        assemble_tetra_nedelec_vector_load_samples_vjp
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_solve_constrained, &
        sparse_direct_solve_constrained_jvp, &
        sparse_direct_solve_constrained_vjp
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none

    private

    public :: solve_tetra_nedelec_sampled_state
    public :: solve_tetra_nedelec_sampled_state_jvp
    public :: solve_tetra_nedelec_sampled_state_vjp

contains

    subroutine solve_tetra_nedelec_sampled_state( &
            vertices, tetrahedra, order, quadrature_degree, curl_coefficient, &
            mass_coefficient, source_values, constrained, constrained_values, &
            state, status)
        real(dp), intent(in) :: vertices(:, :), source_values(:, :, :)
        integer, intent(in) :: tetrahedra(:, :), order, quadrature_degree
        real(dp), intent(in) :: curl_coefficient, mass_coefficient
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: constrained_values(:)
        real(dp), intent(out) :: state(:)
        integer, intent(out) :: status

        type(csc_t) :: matrix
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: load(:)

        state = 0.0_dp
        status = 1
        if (curl_coefficient == 0.0_dp .and. mass_coefficient == 0.0_dp) return
        call assemble_tetra_nedelec_curl_mass_csc( &
            vertices, tetrahedra, matrix, sparse_status, curl_coefficient, &
            mass_coefficient, order)
        if (sparse_status%code /= 0) return
        call assemble_tetra_nedelec_vector_load_samples( &
            vertices, tetrahedra, order, quadrature_degree, source_values, &
            load, sparse_status)
        if (sparse_status%code /= 0) return
        if (.not. valid_state_shapes( &
            matrix, load, constrained, constrained_values, state)) return
        call sparse_direct_solve_constrained( &
            matrix, load, constrained, constrained_values, state, status)
    end subroutine solve_tetra_nedelec_sampled_state

    subroutine solve_tetra_nedelec_sampled_state_jvp( &
            vertices, tetrahedra, order, quadrature_degree, curl_coefficient, &
            mass_coefficient, source_values, source_gradients, constrained, &
            constrained_values, vertices_dot, curl_coefficient_dot, &
            mass_coefficient_dot, source_parameter_dot, &
            constrained_values_dot, state_dot, status)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        integer, intent(in) :: tetrahedra(:, :), order, quadrature_degree
        real(dp), intent(in) :: curl_coefficient, mass_coefficient
        real(dp), intent(in) :: curl_coefficient_dot, mass_coefficient_dot
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        real(dp), intent(in) :: source_parameter_dot(:, :, :)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: constrained_values(:)
        real(dp), intent(in) :: constrained_values_dot(:)
        real(dp), intent(out) :: state_dot(:)
        integer, intent(out) :: status

        type(csc_t) :: matrix, matrix_dot
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: load(:), load_dot(:)

        state_dot = 0.0_dp
        status = 1
        if (curl_coefficient == 0.0_dp .and. mass_coefficient == 0.0_dp) return
        call assemble_tetra_nedelec_curl_mass_csc( &
            vertices, tetrahedra, matrix, sparse_status, curl_coefficient, &
            mass_coefficient, order)
        if (sparse_status%code /= 0) return
        call assemble_tetra_nedelec_curl_mass_csc_jvp( &
            vertices, tetrahedra, order, curl_coefficient, mass_coefficient, &
            vertices_dot, curl_coefficient_dot, mass_coefficient_dot, &
            matrix_dot, sparse_status)
        if (sparse_status%code /= 0) return
        call assemble_tetra_nedelec_vector_load_samples( &
            vertices, tetrahedra, order, quadrature_degree, source_values, &
            load, sparse_status)
        if (sparse_status%code /= 0) return
        call assemble_tetra_nedelec_vector_load_samples_jvp( &
            vertices, tetrahedra, order, quadrature_degree, source_values, &
            source_gradients, vertices_dot, source_parameter_dot, load_dot, &
            sparse_status)
        if (sparse_status%code /= 0) return
        if (.not. valid_state_shapes( &
            matrix, load, constrained, constrained_values, state_dot)) return
        if (size(constrained_values_dot) /= size(constrained_values)) return
        call sparse_direct_solve_constrained_jvp( &
            matrix, load, constrained, constrained_values, matrix_dot%val, &
            load_dot, constrained_values_dot, state_dot, status)
    end subroutine solve_tetra_nedelec_sampled_state_jvp

    subroutine solve_tetra_nedelec_sampled_state_vjp( &
            vertices, tetrahedra, order, quadrature_degree, curl_coefficient, &
            mass_coefficient, source_values, source_gradients, constrained, &
            constrained_values, state, state_bar, vertices_bar, &
            curl_coefficient_bar, mass_coefficient_bar, source_values_bar, &
            constrained_values_bar, status)
        real(dp), intent(in) :: vertices(:, :), source_values(:, :, :)
        integer, intent(in) :: tetrahedra(:, :), order, quadrature_degree
        real(dp), intent(in) :: curl_coefficient, mass_coefficient
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: constrained_values(:), state(:), state_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :)
        real(dp), intent(out) :: curl_coefficient_bar, mass_coefficient_bar
        real(dp), intent(out) :: source_values_bar(:, :, :)
        real(dp), intent(out) :: constrained_values_bar(:)
        integer, intent(out) :: status

        type(csc_t) :: matrix
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: load(:), load_bar(:), matrix_values_bar(:)
        real(dp), allocatable :: load_vertices_bar(:, :)
        real(dp), allocatable :: matrix_vertices_bar(:, :)

        vertices_bar = 0.0_dp
        curl_coefficient_bar = 0.0_dp
        mass_coefficient_bar = 0.0_dp
        source_values_bar = 0.0_dp
        constrained_values_bar = 0.0_dp
        status = 1
        if (curl_coefficient == 0.0_dp .and. mass_coefficient == 0.0_dp) return
        if (any(shape(vertices_bar) /= shape(vertices))) return
        call assemble_tetra_nedelec_curl_mass_csc( &
            vertices, tetrahedra, matrix, sparse_status, curl_coefficient, &
            mass_coefficient, order)
        if (sparse_status%code /= 0) return
        call assemble_tetra_nedelec_vector_load_samples( &
            vertices, tetrahedra, order, quadrature_degree, source_values, &
            load, sparse_status)
        if (sparse_status%code /= 0) return
        if (.not. valid_state_shapes( &
            matrix, load, constrained, constrained_values, state)) return
        if (size(state_bar) /= size(state)) return
        if (size(constrained_values_bar) /= size(constrained_values)) return
        allocate(matrix_values_bar(matrix%nnz), load_bar(size(load)))
        call sparse_direct_solve_constrained_vjp( &
            matrix, load, constrained, constrained_values, state, state_bar, &
            matrix_values_bar, load_bar, constrained_values_bar, status)
        if (status /= 0) return
        allocate(matrix_vertices_bar, mold=vertices)
        allocate(load_vertices_bar, mold=vertices)
        call assemble_tetra_nedelec_curl_mass_csc_vjp( &
            vertices, tetrahedra, order, curl_coefficient, mass_coefficient, &
            matrix_values_bar, matrix_vertices_bar, curl_coefficient_bar, &
            mass_coefficient_bar, sparse_status)
        if (sparse_status%code /= 0) return
        call assemble_tetra_nedelec_vector_load_samples_vjp( &
            vertices, tetrahedra, order, quadrature_degree, source_values, &
            source_gradients, load_bar, load_vertices_bar, source_values_bar, &
            sparse_status)
        if (sparse_status%code /= 0) return
        vertices_bar = matrix_vertices_bar + load_vertices_bar
        status = 0
    end subroutine solve_tetra_nedelec_sampled_state_vjp

    pure logical function valid_state_shapes( &
            matrix, load, constrained, constrained_values, state) result(valid)
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: load(:), constrained_values(:), state(:)
        logical, intent(in) :: constrained(:)

        valid = size(load) == matrix%nrow
        if (.not. valid) return
        valid = size(constrained) == matrix%nrow
        if (.not. valid) return
        valid = size(constrained_values) == matrix%nrow
        if (.not. valid) return
        valid = size(state) == matrix%nrow
    end function valid_state_shapes

end module fortfem_tetra_nedelec_sampled_state_3d
