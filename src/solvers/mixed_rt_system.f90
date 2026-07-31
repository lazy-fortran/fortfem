module fortfem_mixed_rt_system
    !! Differentiable algebraic RT--DG saddle system.
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_factor_csc, &
        sparse_direct_factor_t, sparse_direct_factor_transpose_csc, &
        sparse_direct_free, sparse_direct_solve_factored, &
        sparse_direct_solve_factored_jvp, sparse_direct_solve_factored_vjp
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    private

    public :: solve_mixed_rt_system
    public :: solve_mixed_rt_system_jvp
    public :: solve_mixed_rt_system_vjp

contains

    subroutine solve_mixed_rt_system( &
            flux_mass, divergence, load, flux, pressure, status)
        type(csc_t), intent(in) :: flux_mass, divergence
        real(dp), intent(in) :: load(:)
        real(dp), intent(out) :: flux(:), pressure(:)
        integer, intent(out) :: status

        type(csc_t) :: system
        type(sparse_direct_factor_t) :: factor
        real(dp), allocatable :: right_hand_side(:), state(:)
        integer :: flux_count, pressure_count, system_size

        flux = 0.0_dp
        pressure = 0.0_dp
        call validate_inputs( &
            flux_mass, divergence, load, flux, pressure, status)
        if (status /= 0) return
        flux_count = flux_mass%nrow
        pressure_count = divergence%nrow
        system_size = flux_count + pressure_count
        call build_system(flux_mass, divergence, system, status)
        if (status /= 0) return
        allocate(right_hand_side(system_size), state(system_size))
        right_hand_side = 0.0_dp
        right_hand_side(flux_count + 1:) = load
        call sparse_direct_factor_csc( &
            factor, system_size, system%col_ptr, system%row_idx, system%val, &
            status)
        if (status /= 0) return
        call sparse_direct_solve_factored( &
            factor, right_hand_side, state, status)
        call sparse_direct_free(factor)
        if (status /= 0) return
        flux = state(:flux_count)
        pressure = state(flux_count + 1:)
    end subroutine solve_mixed_rt_system

    subroutine solve_mixed_rt_system_jvp( &
            flux_mass, divergence, load, flux_mass_dot, divergence_dot, &
            load_dot, flux_dot, pressure_dot, status)
        type(csc_t), intent(in) :: flux_mass, divergence
        type(csc_t), intent(in) :: flux_mass_dot, divergence_dot
        real(dp), intent(in) :: load(:), load_dot(:)
        real(dp), intent(out) :: flux_dot(:), pressure_dot(:)
        integer, intent(out) :: status

        type(csc_t) :: system, system_dot
        type(sparse_direct_factor_t) :: factor
        real(dp), allocatable :: right_hand_side(:), right_hand_side_dot(:)
        real(dp), allocatable :: state(:), state_dot(:)
        integer :: flux_count, pressure_count, system_size

        flux_dot = 0.0_dp
        pressure_dot = 0.0_dp
        call validate_inputs( &
            flux_mass, divergence, load, flux_dot, pressure_dot, status)
        if (status /= 0) return
        status = 1
        if (.not. same_pattern(flux_mass, flux_mass_dot)) return
        if (.not. same_pattern(divergence, divergence_dot)) return
        if (size(load_dot) /= size(load)) return
        flux_count = flux_mass%nrow
        pressure_count = divergence%nrow
        system_size = flux_count + pressure_count
        call build_system(flux_mass, divergence, system, status)
        if (status /= 0) return
        call build_system(flux_mass_dot, divergence_dot, system_dot, status)
        if (status /= 0) return
        if (.not. same_pattern(system, system_dot)) then
            status = 1
            return
        end if
        allocate(right_hand_side(system_size), state(system_size))
        allocate(right_hand_side_dot(system_size), state_dot(system_size))
        right_hand_side = 0.0_dp
        right_hand_side_dot = 0.0_dp
        right_hand_side(flux_count + 1:) = load
        right_hand_side_dot(flux_count + 1:) = load_dot
        call sparse_direct_factor_csc( &
            factor, system_size, system%col_ptr, system%row_idx, system%val, &
            status)
        if (status /= 0) return
        call sparse_direct_solve_factored( &
            factor, right_hand_side, state, status)
        if (status == 0) then
            call sparse_direct_solve_factored_jvp( &
                factor, system_size, system%col_ptr, system%row_idx, &
                system_dot%val, state, right_hand_side_dot, state_dot, status)
        end if
        call sparse_direct_free(factor)
        if (status /= 0) return
        flux_dot = state_dot(:flux_count)
        pressure_dot = state_dot(flux_count + 1:)
    end subroutine solve_mixed_rt_system_jvp

    subroutine solve_mixed_rt_system_vjp( &
            flux_mass, divergence, load, flux, pressure, flux_bar, &
            pressure_bar, flux_mass_values_bar, divergence_values_bar, &
            load_bar, status)
        type(csc_t), intent(in) :: flux_mass, divergence
        real(dp), intent(in) :: load(:), flux(:), pressure(:)
        real(dp), intent(in) :: flux_bar(:), pressure_bar(:)
        real(dp), intent(out) :: flux_mass_values_bar(:)
        real(dp), intent(out) :: divergence_values_bar(:), load_bar(:)
        integer, intent(out) :: status

        type(csc_t) :: system
        type(sparse_direct_factor_t) :: transpose_factor
        real(dp), allocatable :: right_hand_side_bar(:)
        real(dp), allocatable :: state(:), state_bar(:), system_values_bar(:)
        integer :: flux_count, pressure_count, system_size

        flux_mass_values_bar = 0.0_dp
        divergence_values_bar = 0.0_dp
        load_bar = 0.0_dp
        call validate_inputs( &
            flux_mass, divergence, load, flux, pressure, status)
        if (status /= 0) return
        status = 1
        if (size(flux_bar) /= size(flux) .or. &
            size(pressure_bar) /= size(pressure)) return
        if (size(flux_mass_values_bar) /= flux_mass%nnz .or. &
            size(divergence_values_bar) /= divergence%nnz .or. &
            size(load_bar) /= size(load)) return
        flux_count = flux_mass%nrow
        pressure_count = divergence%nrow
        system_size = flux_count + pressure_count
        call build_system(flux_mass, divergence, system, status)
        if (status /= 0) return
        allocate(state(system_size), state_bar(system_size))
        allocate(right_hand_side_bar(system_size))
        allocate(system_values_bar(system%nnz))
        state(:flux_count) = flux
        state(flux_count + 1:) = pressure
        state_bar(:flux_count) = flux_bar
        state_bar(flux_count + 1:) = pressure_bar
        call sparse_direct_factor_transpose_csc( &
            transpose_factor, system_size, system%col_ptr, system%row_idx, &
            system%val, status)
        if (status /= 0) return
        call sparse_direct_solve_factored_vjp( &
            transpose_factor, system_size, system%col_ptr, system%row_idx, &
            state, state_bar, right_hand_side_bar, system_values_bar, status)
        call sparse_direct_free(transpose_factor)
        if (status /= 0) return
        load_bar = right_hand_side_bar(flux_count + 1:)
        call gather_block_bar( &
            flux_mass, 0, 0, 1.0_dp, system, system_values_bar, &
            flux_mass_values_bar)
        call gather_block_bar( &
            divergence, flux_count, 0, 1.0_dp, system, system_values_bar, &
            divergence_values_bar)
        call gather_transpose_block_bar( &
            divergence, flux_count, -1.0_dp, system, system_values_bar, &
            divergence_values_bar)
    end subroutine solve_mixed_rt_system_vjp

    subroutine build_system(flux_mass, divergence, system, status)
        type(csc_t), intent(in) :: flux_mass, divergence
        type(csc_t), intent(out) :: system
        integer, intent(out) :: status

        type(fortsparse_status_t) :: sparse_status
        integer, allocatable :: columns(:), rows(:)
        real(dp), allocatable :: values(:)
        integer :: column, entry, flux_count, matrix_entry, system_size

        status = 1
        flux_count = flux_mass%nrow
        system_size = flux_count + divergence%nrow
        allocate(rows(flux_mass%nnz + 2*divergence%nnz))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        call append_block(flux_mass, 0, 0, 1.0_dp)
        call append_block(divergence, flux_count, 0, 1.0_dp)
        call append_transpose_block(divergence, flux_count, -1.0_dp)
        call csc_from_triplet( &
            system_size, system_size, rows, columns, values, system, &
            sparse_status)
        status = sparse_status%code

    contains

        subroutine append_block(block, row_offset, column_offset, scale)
            type(csc_t), intent(in) :: block
            integer, intent(in) :: row_offset, column_offset
            real(dp), intent(in) :: scale

            do column = 1, block%ncol
                do matrix_entry = block%col_ptr(column), &
                        block%col_ptr(column + 1) - 1
                    entry = entry + 1
                    rows(entry) = block%row_idx(matrix_entry) + row_offset
                    columns(entry) = column + column_offset
                    values(entry) = scale*block%val(matrix_entry)
                end do
            end do
        end subroutine append_block

        subroutine append_transpose_block(block, column_offset, scale)
            type(csc_t), intent(in) :: block
            integer, intent(in) :: column_offset
            real(dp), intent(in) :: scale

            do column = 1, block%ncol
                do matrix_entry = block%col_ptr(column), &
                        block%col_ptr(column + 1) - 1
                    entry = entry + 1
                    rows(entry) = column
                    columns(entry) = &
                        block%row_idx(matrix_entry) + column_offset
                    values(entry) = scale*block%val(matrix_entry)
                end do
            end do
        end subroutine append_transpose_block

    end subroutine build_system

    subroutine gather_block_bar( &
            block, row_offset, column_offset, scale, system, &
            system_values_bar, values_bar)
        type(csc_t), intent(in) :: block
        integer, intent(in) :: row_offset, column_offset
        real(dp), intent(in) :: scale
        type(csc_t), intent(in) :: system
        real(dp), intent(in) :: system_values_bar(:)
        real(dp), intent(inout) :: values_bar(:)

        integer :: column, entry

        do column = 1, block%ncol
            do entry = block%col_ptr(column), block%col_ptr(column + 1) - 1
                values_bar(entry) = values_bar(entry) + scale*value_bar_at( &
                    system, system_values_bar, &
                    block%row_idx(entry) + row_offset, column + column_offset)
            end do
        end do
    end subroutine gather_block_bar

    subroutine gather_transpose_block_bar( &
            block, column_offset, scale, system, system_values_bar, values_bar)
        type(csc_t), intent(in) :: block
        integer, intent(in) :: column_offset
        real(dp), intent(in) :: scale
        type(csc_t), intent(in) :: system
        real(dp), intent(in) :: system_values_bar(:)
        real(dp), intent(inout) :: values_bar(:)

        integer :: column, entry

        do column = 1, block%ncol
            do entry = block%col_ptr(column), block%col_ptr(column + 1) - 1
                values_bar(entry) = values_bar(entry) + scale*value_bar_at( &
                    system, system_values_bar, column, &
                    block%row_idx(entry) + column_offset)
            end do
        end do
    end subroutine gather_transpose_block_bar

    pure real(dp) function value_bar_at( &
            system, system_values_bar, row, column) result(value_bar)
        type(csc_t), intent(in) :: system
        real(dp), intent(in) :: system_values_bar(:)
        integer, intent(in) :: row, column

        integer :: entry

        value_bar = 0.0_dp
        do entry = system%col_ptr(column), system%col_ptr(column + 1) - 1
            if (system%row_idx(entry) == row) then
                value_bar = system_values_bar(entry)
                return
            end if
        end do
    end function value_bar_at

    subroutine validate_inputs( &
            flux_mass, divergence, load, flux, pressure, status)
        type(csc_t), intent(in) :: flux_mass, divergence
        real(dp), intent(in) :: load(:), flux(:), pressure(:)
        integer, intent(out) :: status

        status = 1
        if (flux_mass%nrow /= flux_mass%ncol) return
        if (divergence%ncol /= flux_mass%nrow) return
        if (size(load) /= divergence%nrow) return
        if (size(flux) /= flux_mass%nrow) return
        if (size(pressure) /= divergence%nrow) return
        status = 0
    end subroutine validate_inputs

    pure logical function same_pattern(left, right) result(same)
        type(csc_t), intent(in) :: left, right

        same = left%nrow == right%nrow .and. left%ncol == right%ncol
        if (.not. same) return
        same = left%nnz == right%nnz
        if (.not. same) return
        same = all(left%col_ptr == right%col_ptr)
        if (.not. same) return
        same = all(left%row_idx == right%row_idx)
    end function same_pattern

end module fortfem_mixed_rt_system
