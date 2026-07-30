module fortfem_scalar_helmholtz_pml_3d
    !! P1 tetrahedral scalar Helmholtz solver with diagonal complex stretch.
    use fortfem_cartesian_helmholtz_pml, only: &
        cartesian_scalar_helmholtz_pml_coefficients
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3, inv3
    use fortsparse, only: csc_from_triplet, csc_z_t, fortsparse_status_t, &
        sparse_solve_once
    implicit none
    private

    public :: solve_scalar_helmholtz_pml_p1_3d

contains

    subroutine solve_scalar_helmholtz_pml_p1_3d( &
            vertices, tetrahedra, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, solution, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :)
        complex(dp), intent(in) :: stretch(:, :)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: status

        type(csc_z_t) :: matrix, constrained_matrix
        type(fortsparse_status_t) :: sparse_status
        complex(dp), allocatable :: right_hand_side(:), values(:)
        complex(dp) :: element_matrix(4, 4), gradient_coefficient(3)
        complex(dp) :: mass_coefficient
        integer, allocatable :: columns(:), rows(:)
        logical, allocatable :: constrained(:)
        real(dp) :: determinant, gradients(3, 4), inverse_jacobian(3, 3)
        real(dp) :: jacobian(3, 3), local_mass, volume
        integer :: column, constraint, element, entry, inverse_status
        integer :: kept_entry, local_column, local_row, node, row
        integer :: vertex_count

        status = 1
        if (allocated(solution)) deallocate(solution)
        if (size(vertices, 1) /= 3) return
        vertex_count = size(vertices, 2)
        if (vertex_count < 4) return
        if (size(tetrahedra, 1) /= 4) return
        if (size(tetrahedra, 2) < 1) return
        if (size(stretch, 1) /= 3) return
        if (size(stretch, 2) /= size(tetrahedra, 2)) return
        if (size(volume_load) /= vertex_count) return
        if (size(dirichlet_nodes) /= size(dirichlet_values)) return
        if (wave_number <= 0.0_dp) return
        if (any(tetrahedra < 1) .or. any(tetrahedra > vertex_count)) return
        if (any(dirichlet_nodes < 1) .or. &
            any(dirichlet_nodes > vertex_count)) return

        allocate(rows(16*size(tetrahedra, 2)))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        do element = 1, size(tetrahedra, 2)
            jacobian(:, 1) = vertices(:, tetrahedra(2, element)) - &
                vertices(:, tetrahedra(1, element))
            jacobian(:, 2) = vertices(:, tetrahedra(3, element)) - &
                vertices(:, tetrahedra(1, element))
            jacobian(:, 3) = vertices(:, tetrahedra(4, element)) - &
                vertices(:, tetrahedra(1, element))
            determinant = det3(jacobian)
            if (determinant <= tiny(1.0_dp)) return
            call inv3(jacobian, inverse_jacobian, inverse_status)
            if (inverse_status /= 0) return
            gradients(:, 2:4) = transpose(inverse_jacobian)
            gradients(:, 1) = -sum(gradients(:, 2:4), dim=2)
            volume = determinant/6.0_dp
            call cartesian_scalar_helmholtz_pml_coefficients( &
                stretch(:, element), gradient_coefficient, mass_coefficient, &
                inverse_status)
            if (inverse_status /= 0) return
            do local_column = 1, 4
                do local_row = 1, 4
                    local_mass = volume/20.0_dp
                    if (local_row == local_column) local_mass = volume/10.0_dp
                    element_matrix(local_row, local_column) = volume*sum( &
                        gradient_coefficient*gradients(:, local_row)* &
                        gradients(:, local_column)) - &
                        wave_number**2*mass_coefficient*local_mass
                    entry = entry + 1
                    rows(entry) = tetrahedra(local_row, element)
                    columns(entry) = tetrahedra(local_column, element)
                    values(entry) = element_matrix(local_row, local_column)
                end do
            end do
        end do
        call csc_from_triplet( &
            vertex_count, vertex_count, rows, columns, values, matrix, &
            sparse_status)
        if (sparse_status%code /= 0) return

        allocate(right_hand_side(vertex_count), constrained(vertex_count))
        right_hand_side = volume_load
        constrained = .false.
        do constraint = 1, size(dirichlet_nodes)
            node = dirichlet_nodes(constraint)
            if (constrained(node)) return
            constrained(node) = .true.
            do entry = matrix%col_ptr(node), matrix%col_ptr(node + 1) - 1
                row = matrix%row_idx(entry)
                right_hand_side(row) = right_hand_side(row) - &
                    matrix%val(entry)*dirichlet_values(constraint)
            end do
        end do
        deallocate(rows, columns, values)
        allocate(rows(matrix%nnz + size(dirichlet_nodes)))
        allocate(columns(size(rows)), values(size(rows)))
        kept_entry = 0
        do column = 1, matrix%ncol
            if (constrained(column)) cycle
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                if (constrained(row)) cycle
                kept_entry = kept_entry + 1
                rows(kept_entry) = row
                columns(kept_entry) = column
                values(kept_entry) = matrix%val(entry)
            end do
        end do
        do constraint = 1, size(dirichlet_nodes)
            kept_entry = kept_entry + 1
            node = dirichlet_nodes(constraint)
            rows(kept_entry) = node
            columns(kept_entry) = node
            values(kept_entry) = cmplx(1.0_dp, 0.0_dp, dp)
            right_hand_side(node) = dirichlet_values(constraint)
        end do
        call csc_from_triplet( &
            vertex_count, vertex_count, rows(:kept_entry), &
            columns(:kept_entry), values(:kept_entry), constrained_matrix, &
            sparse_status)
        if (sparse_status%code /= 0) return
        allocate(solution(vertex_count))
        call sparse_solve_once( &
            constrained_matrix, right_hand_side, solution, sparse_status)
        if (sparse_status%code /= 0) then
            status = 2
            return
        end if
        status = 0
    end subroutine solve_scalar_helmholtz_pml_p1_3d

end module fortfem_scalar_helmholtz_pml_3d
