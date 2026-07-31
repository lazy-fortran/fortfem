module fortfem_tetra_mixed_poisson_3d
    use fortfem_forms_simple, only: &
        compile_tetra_mixed_form_csc, compile_tetra_rt_form_csc, form_expr_t
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_solve_csc
    use fortfem_tetra_discontinuous_arbitrary_order, only: &
        evaluate_tetra_discontinuous, initialize_tetra_discontinuous, &
        tetra_discontinuous_dof_count, tetra_discontinuous_t
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortnum_linalg, only: det3, det3_jvp, det3_vjp
    use fortsparse, only: &
        csc_from_triplet, csc_t, FORTSPARSE_INTERNAL_ERROR, &
        FORTSPARSE_INVALID_MATRIX, fortsparse_status_t, status_set
    implicit none
    private

    public :: solve_symbolic_tetra_mixed_poisson_rt
    public :: assemble_tetra_dg_source_load_samples
    public :: assemble_tetra_dg_source_load_samples_jvp
    public :: assemble_tetra_dg_source_load_samples_vjp

    abstract interface
        pure function scalar_source_3d(x, y, z) result(value)
            import :: dp
            real(dp), intent(in) :: x, y, z
            real(dp) :: value
        end function scalar_source_3d
    end interface

contains

    subroutine assemble_tetra_dg_source_load_samples( &
            vertices, tetrahedra, degree, quadrature_degree, source_values, &
            load, status)
        real(dp), intent(in) :: vertices(:, :), source_values(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        real(dp), allocatable, intent(out) :: load(:)
        type(fortsparse_status_t), intent(out) :: status

        type(tetra_discontinuous_t) :: basis
        real(dp), allocatable :: basis_values(:), weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, jacobian(3, 3), point(3)
        real(dp) :: tetra_vertices(3, 4)
        integer :: dof, dof_count, local_status, node, point_index, tetrahedron

        call initialize_sampled_dg_load( &
            vertices, tetrahedra, degree, quadrature_degree, source_values, &
            basis, x, y, z, weights, dof_count, status)
        if (status%code /= 0) return
        allocate(load(dof_count*size(tetrahedra, 2)), source=0.0_dp)
        allocate(basis_values(dof_count))
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                tetra_vertices(:, node) = &
                    vertices(:, tetrahedra(node, tetrahedron))
            end do
            call tetra_geometry( &
                tetra_vertices, jacobian, determinant, local_status)
            if (local_status /= 0) return
            do point_index = 1, size(weights)
                point = [x(point_index), y(point_index), z(point_index)]
                call evaluate_tetra_discontinuous( &
                    basis, point(1), point(2), point(3), basis_values, &
                    local_status)
                if (local_status /= 0) return
                do dof = 1, dof_count
                    load((tetrahedron - 1)*dof_count + dof) = &
                        load((tetrahedron - 1)*dof_count + dof) + &
                        determinant*weights(point_index)*basis_values(dof)* &
                        source_values(point_index, tetrahedron)
                end do
            end do
        end do
        call status_set(status, 0, "")
    end subroutine assemble_tetra_dg_source_load_samples

    subroutine assemble_tetra_dg_source_load_samples_jvp( &
            vertices, tetrahedra, degree, quadrature_degree, source_values, &
            source_gradients, vertices_dot, source_parameter_dot, load_dot, &
            status)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        real(dp), intent(in) :: source_values(:, :), source_gradients(:, :, :)
        real(dp), intent(in) :: source_parameter_dot(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        real(dp), allocatable, intent(out) :: load_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        type(tetra_discontinuous_t) :: basis
        real(dp), allocatable :: basis_values(:), weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, determinant_dot, jacobian(3, 3)
        real(dp) :: jacobian_dot(3, 3), point(3), point_dot(3), source_dot
        real(dp) :: tetra_vertices(3, 4), tetra_vertices_dot(3, 4)
        integer :: dof, dof_count, local_status, node, point_index, tetrahedron

        call initialize_sampled_dg_load( &
            vertices, tetrahedra, degree, quadrature_degree, source_values, &
            basis, x, y, z, weights, dof_count, status)
        if (status%code /= 0) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (size(source_gradients, 1) /= 3 .or. &
            size(source_gradients, 2) /= size(weights) .or. &
            size(source_gradients, 3) /= size(tetrahedra, 2)) return
        if (any(shape(source_parameter_dot) /= shape(source_values))) return
        allocate(load_dot(dof_count*size(tetrahedra, 2)), source=0.0_dp)
        allocate(basis_values(dof_count))
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                tetra_vertices(:, node) = &
                    vertices(:, tetrahedra(node, tetrahedron))
                tetra_vertices_dot(:, node) = &
                    vertices_dot(:, tetrahedra(node, tetrahedron))
            end do
            call tetra_geometry( &
                tetra_vertices, jacobian, determinant, local_status)
            if (local_status /= 0) return
            call tetra_jacobian(tetra_vertices_dot, jacobian_dot)
            call det3_jvp(jacobian, jacobian_dot, determinant_dot)
            do point_index = 1, size(weights)
                point = [x(point_index), y(point_index), z(point_index)]
                call evaluate_tetra_discontinuous( &
                    basis, point(1), point(2), point(3), basis_values, &
                    local_status)
                if (local_status /= 0) return
                point_dot = tetra_vertices_dot(:, 1) + &
                    matmul(jacobian_dot, point)
                source_dot = source_parameter_dot(point_index, tetrahedron) + &
                    dot_product( &
                    source_gradients(:, point_index, tetrahedron), point_dot)
                do dof = 1, dof_count
                    load_dot((tetrahedron - 1)*dof_count + dof) = &
                        load_dot((tetrahedron - 1)*dof_count + dof) + &
                        weights(point_index)*basis_values(dof)*( &
                        determinant_dot*source_values( &
                        point_index, tetrahedron) + determinant*source_dot)
                end do
            end do
        end do
        call status_set(status, 0, "")
    end subroutine assemble_tetra_dg_source_load_samples_jvp

    subroutine assemble_tetra_dg_source_load_samples_vjp( &
            vertices, tetrahedra, degree, quadrature_degree, source_values, &
            source_gradients, load_bar, vertices_bar, source_values_bar, &
            status)
        real(dp), intent(in) :: vertices(:, :), source_values(:, :)
        real(dp), intent(in) :: source_gradients(:, :, :), load_bar(:)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        real(dp), intent(out) :: vertices_bar(:, :), source_values_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        type(tetra_discontinuous_t) :: basis
        real(dp), allocatable :: basis_values(:), weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, determinant_bar, determinant_jacobian_bar(3, 3)
        real(dp) :: jacobian(3, 3), jacobian_bar(3, 3), point(3), point_bar(3)
        real(dp) :: sample_bar, seed, tetra_vertices(3, 4)
        real(dp) :: tetra_vertices_bar(3, 4)
        integer :: dof, dof_count, local_status, node, point_index, tetrahedron

        vertices_bar = 0.0_dp
        source_values_bar = 0.0_dp
        call initialize_sampled_dg_load( &
            vertices, tetrahedra, degree, quadrature_degree, source_values, &
            basis, x, y, z, weights, dof_count, status)
        if (status%code /= 0) return
        if (any(shape(vertices_bar) /= shape(vertices))) return
        if (any(shape(source_values_bar) /= shape(source_values))) return
        if (size(source_gradients, 1) /= 3 .or. &
            size(source_gradients, 2) /= size(weights) .or. &
            size(source_gradients, 3) /= size(tetrahedra, 2)) return
        if (size(load_bar) /= dof_count*size(tetrahedra, 2)) return
        allocate(basis_values(dof_count))
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                tetra_vertices(:, node) = &
                    vertices(:, tetrahedra(node, tetrahedron))
            end do
            call tetra_geometry( &
                tetra_vertices, jacobian, determinant, local_status)
            if (local_status /= 0) return
            tetra_vertices_bar = 0.0_dp
            jacobian_bar = 0.0_dp
            determinant_bar = 0.0_dp
            do point_index = 1, size(weights)
                point = [x(point_index), y(point_index), z(point_index)]
                call evaluate_tetra_discontinuous( &
                    basis, point(1), point(2), point(3), basis_values, &
                    local_status)
                if (local_status /= 0) return
                seed = dot_product( &
                    load_bar((tetrahedron - 1)*dof_count + 1: &
                        tetrahedron*dof_count), &
                    weights(point_index)*basis_values)
                determinant_bar = determinant_bar + &
                    seed*source_values(point_index, tetrahedron)
                sample_bar = seed*determinant
                source_values_bar(point_index, tetrahedron) = sample_bar
                point_bar = &
                    sample_bar*source_gradients(:, point_index, tetrahedron)
                tetra_vertices_bar(:, 1) = &
                    tetra_vertices_bar(:, 1) + point_bar
                jacobian_bar = jacobian_bar + &
                    spread(point_bar, 2, 3)*spread(point, 1, 3)
            end do
            call det3_vjp( &
                jacobian, determinant_bar, determinant_jacobian_bar)
            call tetra_jacobian_vjp( &
                jacobian_bar + determinant_jacobian_bar, tetra_vertices_bar)
            do node = 1, 4
                vertices_bar(:, tetrahedra(node, tetrahedron)) = &
                    vertices_bar(:, tetrahedra(node, tetrahedron)) + &
                    tetra_vertices_bar(:, node)
            end do
        end do
        call status_set(status, 0, "")
    end subroutine assemble_tetra_dg_source_load_samples_vjp

    subroutine solve_symbolic_tetra_mixed_poisson_rt( &
            vertices, tetrahedra, degree, quadrature_degree, flux_mass_form, &
            pressure_to_flux_form, balance_form, source, flux_dofs, &
            pressure_dofs, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        type(form_expr_t), intent(in) :: flux_mass_form
        type(form_expr_t), intent(in) :: pressure_to_flux_form
        type(form_expr_t), intent(in) :: balance_form
        procedure(scalar_source_3d) :: source
        real(dp), allocatable, intent(out) :: flux_dofs(:)
        real(dp), allocatable, intent(out) :: pressure_dofs(:)
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: balance, flux_mass, pressure_to_flux, system
        integer, allocatable :: columns(:), rows(:)
        real(dp), allocatable :: load(:), right_hand_side(:), solution(:)
        real(dp), allocatable :: values(:)
        integer :: column, entry, flux_count, matrix_entry
        integer :: pressure_count, solve_status, system_size

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Symbolic tetra RT-DG mixed Poisson solve failed")
        if (degree < 0 .or. quadrature_degree < 0) return
        call compile_tetra_rt_form_csc( &
            flux_mass_form, vertices, tetrahedra, degree, quadrature_degree, &
            flux_mass, status)
        if (status%code /= 0) return
        call compile_tetra_mixed_form_csc( &
            pressure_to_flux_form, vertices, tetrahedra, degree, &
            quadrature_degree, pressure_to_flux, status)
        if (status%code /= 0) return
        call compile_tetra_mixed_form_csc( &
            balance_form, vertices, tetrahedra, degree, quadrature_degree, &
            balance, status)
        if (status%code /= 0) return
        call assemble_tetra_dg_source_load( &
            vertices, tetrahedra, degree, quadrature_degree, source, load, &
            status)
        if (status%code /= 0) return

        flux_count = flux_mass%nrow
        pressure_count = balance%nrow
        if (flux_mass%ncol /= flux_count) return
        if (balance%ncol /= flux_count) return
        if (pressure_to_flux%nrow /= flux_count .or. &
            pressure_to_flux%ncol /= pressure_count) return
        if (size(load) /= pressure_count) return
        system_size = flux_count + pressure_count
        allocate(rows( &
            flux_mass%nnz + pressure_to_flux%nnz + balance%nnz))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        call append_block(flux_mass, 0, 0)
        call append_block(pressure_to_flux, 0, flux_count)
        call append_block(balance, flux_count, 0)
        call csc_from_triplet( &
            system_size, system_size, rows, columns, values, system, status)
        if (status%code /= 0) return

        allocate(right_hand_side(system_size), solution(system_size))
        right_hand_side = 0.0_dp
        right_hand_side(flux_count + 1:) = load
        call sparse_direct_solve_csc( &
            system_size, system%col_ptr, system%row_idx, system%val, &
            right_hand_side, solution, solve_status)
        if (solve_status /= 0) then
            call status_set( &
                status, FORTSPARSE_INTERNAL_ERROR, &
                "Symbolic tetra RT-DG sparse solve failed")
            return
        end if
        allocate(flux_dofs(flux_count), pressure_dofs(pressure_count))
        flux_dofs = solution(:flux_count)
        pressure_dofs = solution(flux_count + 1:)
        call status_set(status, 0, "")

    contains

        subroutine append_block(block, row_offset, column_offset)
            type(csc_t), intent(in) :: block
            integer, intent(in) :: row_offset, column_offset

            do column = 1, block%ncol
                do matrix_entry = block%col_ptr(column), &
                        block%col_ptr(column + 1) - 1
                    entry = entry + 1
                    rows(entry) = block%row_idx(matrix_entry) + row_offset
                    columns(entry) = column + column_offset
                    values(entry) = block%val(matrix_entry)
                end do
            end do
        end subroutine append_block

    end subroutine solve_symbolic_tetra_mixed_poisson_rt

    subroutine assemble_tetra_dg_source_load( &
            vertices, tetrahedra, degree, quadrature_degree, source, load, &
            status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        procedure(scalar_source_3d) :: source
        real(dp), allocatable, intent(out) :: load(:)
        type(fortsparse_status_t), intent(out) :: status

        type(tetra_discontinuous_t) :: basis
        real(dp), allocatable :: basis_values(:), weights(:)
        real(dp), allocatable :: x(:), y(:), z(:)
        real(dp) :: determinant, jacobian(3, 3), physical_point(3)
        real(dp) :: point(3), tetra_vertices(3, 4)
        integer :: dof, dof_count, local_status, node, point_index
        integer :: tetrahedron

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Tetrahedral DG source load failed")
        if (size(vertices, 1) /= 3 .or. size(tetrahedra, 1) /= 4) return
        if (degree < 0 .or. quadrature_degree < 0) return
        call initialize_tetra_discontinuous(degree, basis, local_status)
        if (local_status /= 0) return
        call tetra_duffy_quadrature( &
            quadrature_degree, x, y, z, weights, local_status)
        if (local_status /= 0) return
        dof_count = tetra_discontinuous_dof_count(basis)
        allocate(load(dof_count*size(tetrahedra, 2)))
        allocate(basis_values(dof_count))
        load = 0.0_dp
        do tetrahedron = 1, size(tetrahedra, 2)
            do node = 1, 4
                tetra_vertices(:, node) = &
                    vertices(:, tetrahedra(node, tetrahedron))
            end do
            jacobian(:, 1) = tetra_vertices(:, 2) - tetra_vertices(:, 1)
            jacobian(:, 2) = tetra_vertices(:, 3) - tetra_vertices(:, 1)
            jacobian(:, 3) = tetra_vertices(:, 4) - tetra_vertices(:, 1)
            determinant = det3(jacobian)
            if (determinant <= tiny(1.0_dp)) return
            do point_index = 1, size(weights)
                point = [x(point_index), y(point_index), z(point_index)]
                call evaluate_tetra_discontinuous( &
                    basis, point(1), point(2), point(3), basis_values, &
                    local_status)
                if (local_status /= 0) return
                physical_point = &
                    tetra_vertices(:, 1) + matmul(jacobian, point)
                do dof = 1, dof_count
                    load((tetrahedron - 1)*dof_count + dof) = &
                        load((tetrahedron - 1)*dof_count + dof) + &
                        determinant*weights(point_index)*basis_values(dof)* &
                        source( &
                        physical_point(1), physical_point(2), &
                        physical_point(3))
                end do
            end do
        end do
        call status_set(status, 0, "")
    end subroutine assemble_tetra_dg_source_load

    subroutine initialize_sampled_dg_load( &
            vertices, tetrahedra, degree, quadrature_degree, source_values, &
            basis, x, y, z, weights, dof_count, status)
        real(dp), intent(in) :: vertices(:, :), source_values(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        type(tetra_discontinuous_t), intent(out) :: basis
        real(dp), allocatable, intent(out) :: x(:), y(:), z(:), weights(:)
        integer, intent(out) :: dof_count
        type(fortsparse_status_t), intent(out) :: status

        integer :: local_status

        dof_count = 0
        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sampled tetrahedral DG source load failed")
        if (size(vertices, 1) /= 3 .or. size(tetrahedra, 1) /= 4) return
        if (degree < 0 .or. quadrature_degree < 0) return
        if (any(tetrahedra < 1) .or. any(tetrahedra > size(vertices, 2))) &
            return
        call initialize_tetra_discontinuous(degree, basis, local_status)
        if (local_status /= 0) return
        call tetra_duffy_quadrature( &
            quadrature_degree, x, y, z, weights, local_status)
        if (local_status /= 0) return
        if (size(source_values, 1) /= size(weights) .or. &
            size(source_values, 2) /= size(tetrahedra, 2)) return
        dof_count = tetra_discontinuous_dof_count(basis)
        call status_set(status, 0, "")
    end subroutine initialize_sampled_dg_load

    subroutine tetra_geometry(vertices, jacobian, determinant, status)
        real(dp), intent(in) :: vertices(3, 4)
        real(dp), intent(out) :: jacobian(3, 3), determinant
        integer, intent(out) :: status

        call tetra_jacobian(vertices, jacobian)
        determinant = det3(jacobian)
        status = 1
        if (determinant <= tiny(1.0_dp)) return
        status = 0
    end subroutine tetra_geometry

    pure subroutine tetra_jacobian(vertices, jacobian)
        real(dp), intent(in) :: vertices(3, 4)
        real(dp), intent(out) :: jacobian(3, 3)

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
    end subroutine tetra_jacobian

    pure subroutine tetra_jacobian_vjp(jacobian_bar, vertices_bar)
        real(dp), intent(in) :: jacobian_bar(3, 3)
        real(dp), intent(inout) :: vertices_bar(3, 4)

        vertices_bar(:, 1) = vertices_bar(:, 1) - &
            sum(jacobian_bar, dim=2)
        vertices_bar(:, 2) = vertices_bar(:, 2) + jacobian_bar(:, 1)
        vertices_bar(:, 3) = vertices_bar(:, 3) + jacobian_bar(:, 2)
        vertices_bar(:, 4) = vertices_bar(:, 4) + jacobian_bar(:, 3)
    end subroutine tetra_jacobian_vjp

end module fortfem_tetra_mixed_poisson_3d
