module fortfem_block_graph_csc
    !! Sparse CSC materialization at the retained-factor boundary.
    !!
    !! A packed block graph remains the preferred matrix-free representation.
    !! These adapters materialize only sparse triplets when a caller explicitly
    !! requests a FortSparse CSC matrix for a retained direct factorization.
    !! Duplicate graph edges are passed through to FortSparse and summed there;
    !! no dense global matrix is formed.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_t, csc_z_t, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, fortsparse_status_t, status_set
    implicit none
    private

    public :: assemble_block_graph_csc
    public :: assemble_block_graph_csc_real
    public :: assemble_block_graph_csc_complex

    interface assemble_block_graph_csc
        module procedure assemble_block_graph_csc_real
        module procedure assemble_block_graph_csc_complex
    end interface assemble_block_graph_csc

contains

    subroutine assemble_block_graph_csc_real( &
            field_sizes, edge_rows, edge_columns, block_offsets, block_values, &
            matrix, status)
        integer, intent(in) :: field_sizes(:), edge_rows(:), edge_columns(:)
        integer, intent(in) :: block_offsets(:)
        real(dp), intent(in) :: block_values(:)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: rows(:), columns(:)
        real(dp), allocatable :: values(:)
        integer :: edge, row, column, row_offset, column_offset, value_index
        integer :: entry, total_state, total_entries

        call reset_real_matrix(matrix)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "packed block graph CSC received incompatible arrays")
        if (.not. valid_graph_inputs( &
            field_sizes, edge_rows, edge_columns, block_offsets, size(block_values))) return
        if (.not. all(ieee_is_finite(block_values))) return

        total_state = sum(field_sizes)
        total_entries = graph_entry_count(field_sizes, edge_rows, edge_columns)
        allocate(rows(total_entries), columns(total_entries), values(total_entries))
        entry = 0
        do edge = 1, size(edge_rows)
            row_offset = field_offset(field_sizes, edge_rows(edge))
            column_offset = field_offset(field_sizes, edge_columns(edge))
            do column = 1, field_sizes(edge_columns(edge))
                do row = 1, field_sizes(edge_rows(edge))
                    entry = entry + 1
                    value_index = block_offsets(edge) + row - 1 + &
                        (column - 1)*field_sizes(edge_rows(edge))
                    rows(entry) = row_offset + row - 1
                    columns(entry) = column_offset + column - 1
                    values(entry) = block_values(value_index)
                end do
            end do
        end do
        call csc_from_triplet(total_state, total_state, rows, columns, values, matrix, status)
        if (status%code == FORTSPARSE_OK) return
        call reset_real_matrix(matrix)
    end subroutine assemble_block_graph_csc_real

    subroutine assemble_block_graph_csc_complex( &
            field_sizes, edge_rows, edge_columns, block_offsets, block_values, &
            matrix, status)
        integer, intent(in) :: field_sizes(:), edge_rows(:), edge_columns(:)
        integer, intent(in) :: block_offsets(:)
        complex(dp), intent(in) :: block_values(:)
        type(csc_z_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: rows(:), columns(:)
        complex(dp), allocatable :: values(:)
        integer :: edge, row, column, row_offset, column_offset, value_index
        integer :: entry, total_state, total_entries

        call reset_complex_matrix(matrix)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "complex packed block graph CSC received incompatible arrays")
        if (.not. valid_graph_inputs( &
            field_sizes, edge_rows, edge_columns, block_offsets, size(block_values))) return
        if (.not. finite_complex(block_values)) return

        total_state = sum(field_sizes)
        total_entries = graph_entry_count(field_sizes, edge_rows, edge_columns)
        allocate(rows(total_entries), columns(total_entries), values(total_entries))
        entry = 0
        do edge = 1, size(edge_rows)
            row_offset = field_offset(field_sizes, edge_rows(edge))
            column_offset = field_offset(field_sizes, edge_columns(edge))
            do column = 1, field_sizes(edge_columns(edge))
                do row = 1, field_sizes(edge_rows(edge))
                    entry = entry + 1
                    value_index = block_offsets(edge) + row - 1 + &
                        (column - 1)*field_sizes(edge_rows(edge))
                    rows(entry) = row_offset + row - 1
                    columns(entry) = column_offset + column - 1
                    values(entry) = block_values(value_index)
                end do
            end do
        end do
        call csc_from_triplet(total_state, total_state, rows, columns, values, matrix, status)
        if (status%code == FORTSPARSE_OK) return
        call reset_complex_matrix(matrix)
    end subroutine assemble_block_graph_csc_complex

    logical function valid_graph_inputs( &
            field_sizes, edge_rows, edge_columns, block_offsets, value_count) result(valid)
        integer, intent(in) :: field_sizes(:), edge_rows(:), edge_columns(:)
        integer, intent(in) :: block_offsets(:), value_count
        integer :: edge, expected_values

        valid = .false.
        if (size(field_sizes) < 1 .or. size(edge_rows) < 1 .or. &
            size(edge_columns) /= size(edge_rows) .or. &
            size(block_offsets) /= size(edge_rows) + 1) return
        if (any(field_sizes < 1) .or. value_count < 1) return
        do edge = 1, size(edge_rows)
            if (edge_rows(edge) < 1 .or. edge_rows(edge) > size(field_sizes) .or. &
                edge_columns(edge) < 1 .or. edge_columns(edge) > size(field_sizes)) return
            if (block_offsets(edge) < 1 .or. block_offsets(edge + 1) <= block_offsets(edge)) return
            expected_values = field_sizes(edge_rows(edge))*field_sizes(edge_columns(edge))
            if (block_offsets(edge + 1) - block_offsets(edge) /= expected_values) return
        end do
        if (block_offsets(1) /= 1 .or. block_offsets(size(block_offsets)) /= value_count + 1) return
        valid = .true.
    end function valid_graph_inputs

    integer function graph_entry_count(field_sizes, edge_rows, edge_columns) result(count)
        integer, intent(in) :: field_sizes(:), edge_rows(:), edge_columns(:)
        integer :: edge

        count = 0
        do edge = 1, size(edge_rows)
            count = count + field_sizes(edge_rows(edge))*field_sizes(edge_columns(edge))
        end do
    end function graph_entry_count

    integer function field_offset(field_sizes, field_index) result(offset)
        integer, intent(in) :: field_sizes(:), field_index
        integer :: index

        offset = 1
        do index = 1, field_index - 1
            offset = offset + field_sizes(index)
        end do
    end function field_offset

    logical function finite_complex(values) result(valid)
        complex(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex

    subroutine reset_real_matrix(matrix)
        type(csc_t), intent(out) :: matrix

        matrix%nrow = 0
        matrix%ncol = 0
        matrix%nnz = 0
    end subroutine reset_real_matrix

    subroutine reset_complex_matrix(matrix)
        type(csc_z_t), intent(out) :: matrix

        matrix%nrow = 0
        matrix%ncol = 0
        matrix%nnz = 0
    end subroutine reset_complex_matrix

end module fortfem_block_graph_csc
