module fortfem_triangle_full_vector_sampled_state_2d
    !! Differentiable constrained full-vector triangle states with sampled forcing.
    use fortfem_assembly_full_vector_arbitrary_order_2d, only: &
        assemble_triangle_bdm_div_mass_csc, &
        assemble_triangle_bdm_div_mass_csc_jvp, &
        assemble_triangle_bdm_div_mass_csc_vjp, &
        assemble_triangle_bdm_vector_load_samples, &
        assemble_triangle_bdm_vector_load_samples_jvp, &
        assemble_triangle_bdm_vector_load_samples_vjp, &
        assemble_triangle_nedelec_second_curl_mass_csc, &
        assemble_triangle_nedelec_second_curl_mass_csc_jvp, &
        assemble_triangle_nedelec_second_curl_mass_csc_vjp, &
        assemble_triangle_nedelec_second_vector_load_samples, &
        assemble_triangle_nedelec_second_vector_load_samples_jvp, &
        assemble_triangle_nedelec_second_vector_load_samples_vjp
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_sparse_direct, only: sparse_direct_solve_constrained, &
        sparse_direct_solve_constrained_jvp, &
        sparse_direct_solve_constrained_vjp
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none

    private

    public :: solve_triangle_bdm_sampled_state
    public :: solve_triangle_bdm_sampled_state_jvp
    public :: solve_triangle_bdm_sampled_state_vjp
    public :: solve_triangle_nedelec_second_sampled_state
    public :: solve_triangle_nedelec_second_sampled_state_jvp
    public :: solve_triangle_nedelec_second_sampled_state_vjp

contains

    subroutine solve_triangle_bdm_sampled_state( &
            mesh, degree, quadrature_degree, derivative_coefficient, &
            mass_coefficient, source_values, constrained, constrained_values, &
            state, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(in) :: source_values(:, :, :)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: constrained_values(:)
        real(dp), intent(out) :: state(:)
        integer, intent(out) :: status

        call solve_triangle_full_vector_sampled_state( &
            mesh, degree, quadrature_degree, .true., derivative_coefficient, &
            mass_coefficient, source_values, constrained, constrained_values, &
            state, status)
    end subroutine solve_triangle_bdm_sampled_state

    subroutine solve_triangle_nedelec_second_sampled_state( &
            mesh, degree, quadrature_degree, derivative_coefficient, &
            mass_coefficient, source_values, constrained, constrained_values, &
            state, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(in) :: source_values(:, :, :)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: constrained_values(:)
        real(dp), intent(out) :: state(:)
        integer, intent(out) :: status

        call solve_triangle_full_vector_sampled_state( &
            mesh, degree, quadrature_degree, .false., derivative_coefficient, &
            mass_coefficient, source_values, constrained, constrained_values, &
            state, status)
    end subroutine solve_triangle_nedelec_second_sampled_state

    subroutine solve_triangle_bdm_sampled_state_jvp( &
            mesh, degree, quadrature_degree, derivative_coefficient, &
            mass_coefficient, source_values, source_gradients, constrained, &
            constrained_values, mesh_vertices_dot, derivative_coefficient_dot, &
            mass_coefficient_dot, source_parameter_dot, &
            constrained_values_dot, state_dot, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: constrained_values(:)
        real(dp), intent(in) :: mesh_vertices_dot(:, :)
        real(dp), intent(in) :: derivative_coefficient_dot
        real(dp), intent(in) :: mass_coefficient_dot
        real(dp), intent(in) :: source_parameter_dot(:, :, :)
        real(dp), intent(in) :: constrained_values_dot(:)
        real(dp), intent(out) :: state_dot(:)
        integer, intent(out) :: status

        call solve_triangle_full_vector_sampled_state_jvp( &
            mesh, degree, quadrature_degree, .true., derivative_coefficient, &
            mass_coefficient, source_values, source_gradients, constrained, &
            constrained_values, mesh_vertices_dot, derivative_coefficient_dot, &
            mass_coefficient_dot, source_parameter_dot, &
            constrained_values_dot, state_dot, status)
    end subroutine solve_triangle_bdm_sampled_state_jvp

    subroutine solve_triangle_nedelec_second_sampled_state_jvp( &
            mesh, degree, quadrature_degree, derivative_coefficient, &
            mass_coefficient, source_values, source_gradients, constrained, &
            constrained_values, mesh_vertices_dot, derivative_coefficient_dot, &
            mass_coefficient_dot, source_parameter_dot, &
            constrained_values_dot, state_dot, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: constrained_values(:)
        real(dp), intent(in) :: mesh_vertices_dot(:, :)
        real(dp), intent(in) :: derivative_coefficient_dot
        real(dp), intent(in) :: mass_coefficient_dot
        real(dp), intent(in) :: source_parameter_dot(:, :, :)
        real(dp), intent(in) :: constrained_values_dot(:)
        real(dp), intent(out) :: state_dot(:)
        integer, intent(out) :: status

        call solve_triangle_full_vector_sampled_state_jvp( &
            mesh, degree, quadrature_degree, .false., derivative_coefficient, &
            mass_coefficient, source_values, source_gradients, constrained, &
            constrained_values, mesh_vertices_dot, derivative_coefficient_dot, &
            mass_coefficient_dot, source_parameter_dot, &
            constrained_values_dot, state_dot, status)
    end subroutine solve_triangle_nedelec_second_sampled_state_jvp

    subroutine solve_triangle_bdm_sampled_state_vjp( &
            mesh, degree, quadrature_degree, derivative_coefficient, &
            mass_coefficient, source_values, source_gradients, constrained, &
            constrained_values, state, state_bar, mesh_vertices_bar, &
            derivative_coefficient_bar, mass_coefficient_bar, &
            source_values_bar, constrained_values_bar, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: constrained_values(:), state(:), state_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        real(dp), intent(out) :: derivative_coefficient_bar
        real(dp), intent(out) :: mass_coefficient_bar
        real(dp), intent(out) :: source_values_bar(:, :, :)
        real(dp), intent(out) :: constrained_values_bar(:)
        integer, intent(out) :: status

        call solve_triangle_full_vector_sampled_state_vjp( &
            mesh, degree, quadrature_degree, .true., derivative_coefficient, &
            mass_coefficient, source_values, source_gradients, constrained, &
            constrained_values, state, state_bar, mesh_vertices_bar, &
            derivative_coefficient_bar, mass_coefficient_bar, &
            source_values_bar, constrained_values_bar, status)
    end subroutine solve_triangle_bdm_sampled_state_vjp

    subroutine solve_triangle_nedelec_second_sampled_state_vjp( &
            mesh, degree, quadrature_degree, derivative_coefficient, &
            mass_coefficient, source_values, source_gradients, constrained, &
            constrained_values, state, state_bar, mesh_vertices_bar, &
            derivative_coefficient_bar, mass_coefficient_bar, &
            source_values_bar, constrained_values_bar, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: constrained_values(:), state(:), state_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        real(dp), intent(out) :: derivative_coefficient_bar
        real(dp), intent(out) :: mass_coefficient_bar
        real(dp), intent(out) :: source_values_bar(:, :, :)
        real(dp), intent(out) :: constrained_values_bar(:)
        integer, intent(out) :: status

        call solve_triangle_full_vector_sampled_state_vjp( &
            mesh, degree, quadrature_degree, .false., derivative_coefficient, &
            mass_coefficient, source_values, source_gradients, constrained, &
            constrained_values, state, state_bar, mesh_vertices_bar, &
            derivative_coefficient_bar, mass_coefficient_bar, &
            source_values_bar, constrained_values_bar, status)
    end subroutine solve_triangle_nedelec_second_sampled_state_vjp

    subroutine solve_triangle_full_vector_sampled_state( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, source_values, &
            constrained, constrained_values, state, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(in) :: source_values(:, :, :)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: constrained_values(:)
        real(dp), intent(out) :: state(:)
        integer, intent(out) :: status

        type(csc_t) :: matrix
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: load(:)

        state = 0.0_dp
        status = 1
        call assemble_matrix( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, matrix, sparse_status)
        if (sparse_status%code /= 0) return
        call assemble_load( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
            load, sparse_status)
        if (sparse_status%code /= 0) return
        if (.not. valid_state_shapes( &
            matrix, load, constrained, constrained_values, state)) return
        call sparse_direct_solve_constrained( &
            matrix, load, constrained, constrained_values, state, status)
    end subroutine solve_triangle_full_vector_sampled_state

    subroutine solve_triangle_full_vector_sampled_state_jvp( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, source_values, &
            source_gradients, constrained, constrained_values, &
            mesh_vertices_dot, derivative_coefficient_dot, &
            mass_coefficient_dot, source_parameter_dot, &
            constrained_values_dot, state_dot, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: constrained_values(:)
        real(dp), intent(in) :: mesh_vertices_dot(:, :)
        real(dp), intent(in) :: derivative_coefficient_dot
        real(dp), intent(in) :: mass_coefficient_dot
        real(dp), intent(in) :: source_parameter_dot(:, :, :)
        real(dp), intent(in) :: constrained_values_dot(:)
        real(dp), intent(out) :: state_dot(:)
        integer, intent(out) :: status

        type(csc_t) :: matrix, matrix_dot
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: load(:), load_dot(:)

        state_dot = 0.0_dp
        status = 1
        call assemble_matrix( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, matrix, sparse_status)
        if (sparse_status%code /= 0) return
        call assemble_matrix_jvp( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, mesh_vertices_dot, &
            derivative_coefficient_dot, mass_coefficient_dot, matrix_dot, &
            sparse_status)
        if (sparse_status%code /= 0) return
        call assemble_load( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
            load, sparse_status)
        if (sparse_status%code /= 0) return
        call assemble_load_jvp( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
            source_gradients, mesh_vertices_dot, source_parameter_dot, &
            load_dot, sparse_status)
        if (sparse_status%code /= 0) return
        if (.not. valid_state_shapes( &
            matrix, load, constrained, constrained_values, state_dot)) return
        if (size(constrained_values_dot) /= size(constrained_values)) return
        call sparse_direct_solve_constrained_jvp( &
            matrix, load, constrained, constrained_values, matrix_dot%val, &
            load_dot, constrained_values_dot, state_dot, status)
    end subroutine solve_triangle_full_vector_sampled_state_jvp

    subroutine solve_triangle_full_vector_sampled_state_vjp( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, source_values, &
            source_gradients, constrained, constrained_values, state, &
            state_bar, mesh_vertices_bar, derivative_coefficient_bar, &
            mass_coefficient_bar, source_values_bar, constrained_values_bar, &
            status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: constrained_values(:), state(:), state_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        real(dp), intent(out) :: derivative_coefficient_bar
        real(dp), intent(out) :: mass_coefficient_bar
        real(dp), intent(out) :: source_values_bar(:, :, :)
        real(dp), intent(out) :: constrained_values_bar(:)
        integer, intent(out) :: status

        type(csc_t) :: matrix
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: load(:), load_bar(:), matrix_values_bar(:)
        real(dp), allocatable :: load_vertices_bar(:, :)
        real(dp), allocatable :: matrix_vertices_bar(:, :)

        mesh_vertices_bar = 0.0_dp
        derivative_coefficient_bar = 0.0_dp
        mass_coefficient_bar = 0.0_dp
        source_values_bar = 0.0_dp
        constrained_values_bar = 0.0_dp
        status = 1
        if (any(shape(mesh_vertices_bar) /= shape(mesh%vertices))) return
        call assemble_matrix( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, matrix, sparse_status)
        if (sparse_status%code /= 0) return
        call assemble_load( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
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
        allocate(matrix_vertices_bar, mold=mesh%vertices)
        allocate(load_vertices_bar, mold=mesh%vertices)
        call assemble_matrix_vjp( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, matrix_values_bar, &
            matrix_vertices_bar, derivative_coefficient_bar, &
            mass_coefficient_bar, sparse_status)
        if (sparse_status%code /= 0) return
        call assemble_load_vjp( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
            source_gradients, load_bar, load_vertices_bar, source_values_bar, &
            sparse_status)
        if (sparse_status%code /= 0) return
        mesh_vertices_bar = matrix_vertices_bar + load_vertices_bar
        status = 0
    end subroutine solve_triangle_full_vector_sampled_state_vjp

    subroutine assemble_matrix( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, matrix, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        if (normal_family) then
            call assemble_triangle_bdm_div_mass_csc( &
                mesh, degree, quadrature_degree, matrix, status, &
                derivative_coefficient, mass_coefficient)
        else
            call assemble_triangle_nedelec_second_curl_mass_csc( &
                mesh, degree, quadrature_degree, matrix, status, &
                derivative_coefficient, mass_coefficient)
        end if
    end subroutine assemble_matrix

    subroutine assemble_matrix_jvp( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, mesh_vertices_dot, &
            derivative_coefficient_dot, mass_coefficient_dot, matrix_dot, &
            status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(in) :: mesh_vertices_dot(:, :)
        real(dp), intent(in) :: derivative_coefficient_dot
        real(dp), intent(in) :: mass_coefficient_dot
        type(csc_t), intent(out) :: matrix_dot
        type(fortsparse_status_t), intent(out) :: status

        if (normal_family) then
            call assemble_triangle_bdm_div_mass_csc_jvp( &
                mesh, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, mesh_vertices_dot, &
                derivative_coefficient_dot, mass_coefficient_dot, matrix_dot, &
                status)
        else
            call assemble_triangle_nedelec_second_curl_mass_csc_jvp( &
                mesh, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, mesh_vertices_dot, &
                derivative_coefficient_dot, mass_coefficient_dot, matrix_dot, &
                status)
        end if
    end subroutine assemble_matrix_jvp

    subroutine assemble_matrix_vjp( &
            mesh, degree, quadrature_degree, normal_family, &
            derivative_coefficient, mass_coefficient, matrix_values_bar, &
            mesh_vertices_bar, derivative_coefficient_bar, &
            mass_coefficient_bar, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: derivative_coefficient, mass_coefficient
        real(dp), intent(in) :: matrix_values_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        real(dp), intent(out) :: derivative_coefficient_bar
        real(dp), intent(out) :: mass_coefficient_bar
        type(fortsparse_status_t), intent(out) :: status

        if (normal_family) then
            call assemble_triangle_bdm_div_mass_csc_vjp( &
                mesh, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, matrix_values_bar, mesh_vertices_bar, &
                derivative_coefficient_bar, mass_coefficient_bar, status)
        else
            call assemble_triangle_nedelec_second_curl_mass_csc_vjp( &
                mesh, degree, quadrature_degree, derivative_coefficient, &
                mass_coefficient, matrix_values_bar, mesh_vertices_bar, &
                derivative_coefficient_bar, mass_coefficient_bar, status)
        end if
    end subroutine assemble_matrix_vjp

    subroutine assemble_load( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
            load, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), allocatable, intent(out) :: load(:)
        type(fortsparse_status_t), intent(out) :: status

        if (normal_family) then
            call assemble_triangle_bdm_vector_load_samples( &
                mesh, degree, quadrature_degree, source_values, load, status)
        else
            call assemble_triangle_nedelec_second_vector_load_samples( &
                mesh, degree, quadrature_degree, source_values, load, status)
        end if
    end subroutine assemble_load

    subroutine assemble_load_jvp( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
            source_gradients, mesh_vertices_dot, source_parameter_dot, &
            load_dot, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :)
        real(dp), intent(in) :: mesh_vertices_dot(:, :)
        real(dp), intent(in) :: source_parameter_dot(:, :, :)
        real(dp), allocatable, intent(out) :: load_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        if (normal_family) then
            call assemble_triangle_bdm_vector_load_samples_jvp( &
                mesh, degree, quadrature_degree, source_values, &
                source_gradients, mesh_vertices_dot, source_parameter_dot, &
                load_dot, status)
        else
            call assemble_triangle_nedelec_second_vector_load_samples_jvp( &
                mesh, degree, quadrature_degree, source_values, &
                source_gradients, mesh_vertices_dot, source_parameter_dot, &
                load_dot, status)
        end if
    end subroutine assemble_load_jvp

    subroutine assemble_load_vjp( &
            mesh, degree, quadrature_degree, normal_family, source_values, &
            source_gradients, load_bar, mesh_vertices_bar, source_values_bar, &
            status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, quadrature_degree
        logical, intent(in) :: normal_family
        real(dp), intent(in) :: source_values(:, :, :)
        real(dp), intent(in) :: source_gradients(:, :, :, :), load_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        real(dp), intent(out) :: source_values_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        if (normal_family) then
            call assemble_triangle_bdm_vector_load_samples_vjp( &
                mesh, degree, quadrature_degree, source_values, &
                source_gradients, load_bar, mesh_vertices_bar, &
                source_values_bar, status)
        else
            call assemble_triangle_nedelec_second_vector_load_samples_vjp( &
                mesh, degree, quadrature_degree, source_values, &
                source_gradients, load_bar, mesh_vertices_bar, &
                source_values_bar, status)
        end if
    end subroutine assemble_load_vjp

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

end module fortfem_triangle_full_vector_sampled_state_2d
