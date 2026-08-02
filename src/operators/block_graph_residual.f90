module fortfem_block_graph_residual
    !! Bounded-memory N-field residual on a packed rectangular block graph.
    !!
    !! A graph edge e maps field edge_columns(e) to field edge_rows(e).  Its
    !! dense local block is stored column-major in
    !! block_values(block_offsets(e):block_offsets(e+1)-1).  The global
    !! residual is the sum of all edge actions minus rhs.  Duplicate edges are
    !! intentional and are accumulated, so callers can retain element, trace,
    !! FEM, BEM, DtN, PML, tensor, or Fourier ownership without materializing a
    !! monolithic matrix.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: assemble_block_graph_residual
    public :: assemble_block_graph_residual_jvp
    public :: assemble_block_graph_residual_vjp

contains

    subroutine assemble_block_graph_residual( &
            field_sizes, edge_rows, edge_columns, block_offsets, block_values, &
            state, rhs, residual, status)
        integer, intent(in) :: field_sizes(:), edge_rows(:), edge_columns(:)
        integer, intent(in) :: block_offsets(:)
        real(dp), intent(in) :: block_values(:), state(:), rhs(:)
        real(dp), intent(out) :: residual(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: edge, row, column, row_offset, column_offset, value_index

        residual = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "packed block graph residual received incompatible arrays")
        if (.not. valid_inputs(field_sizes, edge_rows, edge_columns, block_offsets, &
            block_values, state, rhs, residual)) return
        residual = -rhs
        do edge = 1, size(edge_rows)
            row_offset = field_offset(field_sizes, edge_rows(edge))
            column_offset = field_offset(field_sizes, edge_columns(edge))
            do column = 1, field_sizes(edge_columns(edge))
                do row = 1, field_sizes(edge_rows(edge))
                    value_index = block_offsets(edge) + row - 1 + &
                        (column - 1)*field_sizes(edge_rows(edge))
                    residual(row_offset + row - 1) = residual(row_offset + row - 1) + &
                        block_values(value_index)*state(column_offset + column - 1)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_block_graph_residual

    subroutine assemble_block_graph_residual_jvp( &
            field_sizes, edge_rows, edge_columns, block_offsets, block_values, state, rhs, &
            block_values_dot, state_dot, rhs_dot, residual_dot, status)
        integer, intent(in) :: field_sizes(:), edge_rows(:), edge_columns(:)
        integer, intent(in) :: block_offsets(:)
        real(dp), intent(in) :: block_values(:), state(:), rhs(:)
        real(dp), intent(in) :: block_values_dot(:), state_dot(:), rhs_dot(:)
        real(dp), intent(out) :: residual_dot(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: edge, row, column, row_offset, column_offset, value_index

        residual_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "packed block graph residual JVP received incompatible arrays")
        if (.not. valid_inputs(field_sizes, edge_rows, edge_columns, block_offsets, &
            block_values, state, rhs, residual_dot)) return
        if (size(block_values_dot) /= size(block_values) .or. &
            size(state_dot) /= size(state) .or. size(rhs_dot) /= size(rhs)) return
        if (.not. finite_vectors(block_values_dot, state_dot, rhs_dot)) return
        residual_dot = -rhs_dot
        do edge = 1, size(edge_rows)
            row_offset = field_offset(field_sizes, edge_rows(edge))
            column_offset = field_offset(field_sizes, edge_columns(edge))
            do column = 1, field_sizes(edge_columns(edge))
                do row = 1, field_sizes(edge_rows(edge))
                    value_index = block_offsets(edge) + row - 1 + &
                        (column - 1)*field_sizes(edge_rows(edge))
                    residual_dot(row_offset + row - 1) = residual_dot(row_offset + row - 1) + &
                        block_values_dot(value_index)*state(column_offset + column - 1) + &
                        block_values(value_index)*state_dot(column_offset + column - 1)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_block_graph_residual_jvp

    subroutine assemble_block_graph_residual_vjp( &
            field_sizes, edge_rows, edge_columns, block_offsets, block_values, state, rhs, &
            residual_bar, block_values_bar, state_bar, rhs_bar, status)
        integer, intent(in) :: field_sizes(:), edge_rows(:), edge_columns(:)
        integer, intent(in) :: block_offsets(:)
        real(dp), intent(in) :: block_values(:), state(:), rhs(:), residual_bar(:)
        real(dp), intent(out) :: block_values_bar(:), state_bar(:), rhs_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: edge, row, column, row_offset, column_offset, value_index

        block_values_bar = 0.0_dp
        state_bar = 0.0_dp
        rhs_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "packed block graph residual VJP received incompatible arrays")
        if (.not. valid_inputs(field_sizes, edge_rows, edge_columns, block_offsets, &
            block_values, state, rhs, residual_bar)) return
        if (size(block_values_bar) /= size(block_values) .or. &
            size(state_bar) /= size(state) .or. size(rhs_bar) /= size(rhs)) return
        if (.not. finite_vectors(residual_bar)) return
        rhs_bar = -residual_bar
        do edge = 1, size(edge_rows)
            row_offset = field_offset(field_sizes, edge_rows(edge))
            column_offset = field_offset(field_sizes, edge_columns(edge))
            do column = 1, field_sizes(edge_columns(edge))
                do row = 1, field_sizes(edge_rows(edge))
                    value_index = block_offsets(edge) + row - 1 + &
                        (column - 1)*field_sizes(edge_rows(edge))
                    block_values_bar(value_index) = block_values_bar(value_index) + &
                        residual_bar(row_offset + row - 1)*state(column_offset + column - 1)
                    state_bar(column_offset + column - 1) = state_bar(column_offset + column - 1) + &
                        block_values(value_index)*residual_bar(row_offset + row - 1)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_block_graph_residual_vjp

    logical function valid_inputs( &
            field_sizes, edge_rows, edge_columns, block_offsets, block_values, state, rhs, &
            residual) result(valid)
        integer, intent(in) :: field_sizes(:), edge_rows(:), edge_columns(:)
        integer, intent(in) :: block_offsets(:)
        real(dp), intent(in) :: block_values(:), state(:), rhs(:), residual(:)
        integer :: edge, total_state, expected_values

        valid = .false.
        if (size(field_sizes) < 1 .or. size(edge_rows) < 1 .or. &
            size(edge_columns) /= size(edge_rows) .or. &
            size(block_offsets) /= size(edge_rows) + 1) return
        if (any(field_sizes < 1) .or. size(block_values) < 1) return
        total_state = sum(field_sizes)
        if (size(state) /= total_state .or. size(rhs) /= total_state .or. &
            size(residual) /= total_state) return
        do edge = 1, size(edge_rows)
            if (edge_rows(edge) < 1 .or. edge_rows(edge) > size(field_sizes)) return
            if (edge_columns(edge) < 1 .or. edge_columns(edge) > size(field_sizes)) return
            if (block_offsets(edge) < 1 .or. block_offsets(edge + 1) <= block_offsets(edge)) return
            expected_values = field_sizes(edge_rows(edge))*field_sizes(edge_columns(edge))
            if (block_offsets(edge + 1) - block_offsets(edge) /= expected_values) return
        end do
        if (block_offsets(1) /= 1 .or. block_offsets(size(block_offsets)) /= &
            size(block_values) + 1) return
        if (.not. finite_vectors(block_values, state, rhs)) return
        valid = .true.
    end function valid_inputs

    logical function finite_vectors(values_1, values_2, values_3) result(valid)
        real(dp), intent(in) :: values_1(:)
        real(dp), intent(in), optional :: values_2(:), values_3(:)

        valid = all(ieee_is_finite(values_1))
        if (present(values_2)) valid = valid .and. all(ieee_is_finite(values_2))
        if (present(values_3)) valid = valid .and. all(ieee_is_finite(values_3))
    end function finite_vectors

    integer function field_offset(field_sizes, field_index) result(offset)
        integer, intent(in) :: field_sizes(:), field_index
        integer :: index

        offset = 1
        do index = 1, field_index - 1
            offset = offset + field_sizes(index)
        end do
    end function field_offset

end module fortfem_block_graph_residual
