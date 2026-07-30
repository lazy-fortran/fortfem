module fortfem_tetra_mixed_poisson_3d
    use fortfem_forms_simple, only: &
        compile_tetra_mixed_form_csc, compile_tetra_rt_form_csc, form_expr_t
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_solve_csc
    use fortfem_tetra_discontinuous_arbitrary_order, only: &
        evaluate_tetra_discontinuous, initialize_tetra_discontinuous, &
        tetra_discontinuous_dof_count, tetra_discontinuous_t
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortnum_linalg, only: det3
    use fortsparse, only: &
        csc_from_triplet, csc_t, FORTSPARSE_INTERNAL_ERROR, &
        FORTSPARSE_INVALID_MATRIX, fortsparse_status_t, status_set
    implicit none
    private

    public :: solve_symbolic_tetra_mixed_poisson_rt

    abstract interface
        pure function scalar_source_3d(x, y, z) result(value)
            import :: dp
            real(dp), intent(in) :: x, y, z
            real(dp) :: value
        end function scalar_source_3d
    end interface

contains

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

end module fortfem_tetra_mixed_poisson_3d
