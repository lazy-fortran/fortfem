module fortfem_planar_maxwell_dtn_system
    !! Differentiable discrete curl--curl system with a planar Maxwell DtN map.
    use fortfem_kinds, only: dp
    use fortfem_planar_nedelec_maxwell_dtn, only: &
        pullback_planar_maxwell_dtn_form, &
        pullback_planar_maxwell_dtn_form_jvp, &
        pullback_planar_maxwell_dtn_form_vjp
    use fortfem_sparse_direct, only: sparse_direct_factor_adjoint_csc, &
        sparse_direct_factor_csc, sparse_direct_factor_t, sparse_direct_free, &
        sparse_direct_solve_factored, sparse_direct_solve_factored_jvp, &
        sparse_direct_solve_factored_vjp
    use fortsparse, only: csc_from_triplet, csc_z_t, fortsparse_status_t
    implicit none
    private

    public :: solve_planar_maxwell_dtn_system
    public :: solve_planar_maxwell_dtn_system_jvp
    public :: solve_planar_maxwell_dtn_system_vjp

contains

    subroutine solve_planar_maxwell_dtn_system( &
            volume_matrix, sampling, nx, ny, wave_number, length_x, length_y, &
            load, state, status)
        complex(dp), intent(in) :: volume_matrix(:, :), sampling(:, :)
        integer, intent(in) :: nx, ny
        real(dp), intent(in) :: wave_number, length_x, length_y
        complex(dp), intent(in) :: load(:)
        complex(dp), intent(out) :: state(:)
        integer, intent(out) :: status

        type(csc_z_t) :: matrix_csc
        type(fortsparse_status_t) :: sparse_status
        type(sparse_direct_factor_t) :: factor
        complex(dp), allocatable :: boundary_form(:, :), matrix(:, :)

        state = cmplx(0.0_dp, 0.0_dp, dp)
        call validate_inputs( &
            volume_matrix, sampling, nx, ny, load, state, status)
        if (status /= 0) return
        allocate(boundary_form(size(state), size(state)))
        call pullback_planar_maxwell_dtn_form( &
            sampling, nx, ny, wave_number, length_x, length_y, boundary_form, &
            status)
        if (status /= 0) return
        matrix = volume_matrix + boundary_form
        call dense_to_csc(matrix, matrix_csc, sparse_status)
        if (sparse_status%code /= 0) return
        call sparse_direct_factor_csc( &
            factor, size(state), matrix_csc%col_ptr, matrix_csc%row_idx, &
            matrix_csc%val, status)
        if (status /= 0) return
        call sparse_direct_solve_factored(factor, load, state, status)
        call sparse_direct_free(factor)
    end subroutine solve_planar_maxwell_dtn_system

    subroutine solve_planar_maxwell_dtn_system_jvp( &
            volume_matrix, sampling, nx, ny, wave_number, length_x, length_y, &
            load, volume_matrix_dot, sampling_dot, wave_number_dot, &
            length_x_dot, length_y_dot, load_dot, state_dot, status)
        complex(dp), intent(in) :: volume_matrix(:, :), sampling(:, :)
        integer, intent(in) :: nx, ny
        real(dp), intent(in) :: wave_number, length_x, length_y
        complex(dp), intent(in) :: load(:)
        complex(dp), intent(in) :: volume_matrix_dot(:, :), sampling_dot(:, :)
        real(dp), intent(in) :: wave_number_dot, length_x_dot, length_y_dot
        complex(dp), intent(in) :: load_dot(:)
        complex(dp), intent(out) :: state_dot(:)
        integer, intent(out) :: status

        type(csc_z_t) :: matrix_csc, matrix_dot_csc
        type(fortsparse_status_t) :: sparse_status
        type(sparse_direct_factor_t) :: factor
        complex(dp), allocatable :: boundary_form(:, :), boundary_form_dot(:, :)
        complex(dp), allocatable :: matrix(:, :), matrix_dot(:, :), state(:)

        state_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call validate_inputs( &
            volume_matrix, sampling, nx, ny, load, state_dot, status)
        if (status /= 0) return
        if (any(shape(volume_matrix_dot) /= shape(volume_matrix)) .or. &
            any(shape(sampling_dot) /= shape(sampling)) .or. &
            size(load_dot) /= size(load)) then
            status = 1
            return
        end if
        allocate(boundary_form(size(load), size(load)))
        allocate(boundary_form_dot(size(load), size(load)), state(size(load)))
        call pullback_planar_maxwell_dtn_form( &
            sampling, nx, ny, wave_number, length_x, length_y, boundary_form, &
            status)
        if (status /= 0) return
        call pullback_planar_maxwell_dtn_form_jvp( &
            sampling, nx, ny, wave_number, length_x, length_y, sampling_dot, &
            wave_number_dot, length_x_dot, length_y_dot, boundary_form_dot, &
            status)
        if (status /= 0) return
        matrix = volume_matrix + boundary_form
        matrix_dot = volume_matrix_dot + boundary_form_dot
        call dense_to_csc(matrix, matrix_csc, sparse_status)
        if (sparse_status%code /= 0) return
        call dense_to_csc(matrix_dot, matrix_dot_csc, sparse_status)
        if (sparse_status%code /= 0) return
        call sparse_direct_factor_csc( &
            factor, size(load), matrix_csc%col_ptr, matrix_csc%row_idx, &
            matrix_csc%val, status)
        if (status /= 0) return
        call sparse_direct_solve_factored(factor, load, state, status)
        if (status == 0) then
            call sparse_direct_solve_factored_jvp( &
                factor, size(load), matrix_dot_csc%col_ptr, &
                matrix_dot_csc%row_idx, matrix_dot_csc%val, state, load_dot, &
                state_dot, status)
        end if
        call sparse_direct_free(factor)
    end subroutine solve_planar_maxwell_dtn_system_jvp

    subroutine solve_planar_maxwell_dtn_system_vjp( &
            volume_matrix, sampling, nx, ny, wave_number, length_x, length_y, &
            load, state, state_bar, volume_matrix_bar, sampling_bar, &
            wave_number_bar, length_x_bar, length_y_bar, load_bar, status)
        complex(dp), intent(in) :: volume_matrix(:, :), sampling(:, :)
        integer, intent(in) :: nx, ny
        real(dp), intent(in) :: wave_number, length_x, length_y
        complex(dp), intent(in) :: load(:), state(:), state_bar(:)
        complex(dp), intent(out) :: volume_matrix_bar(:, :), sampling_bar(:, :)
        real(dp), intent(out) :: wave_number_bar, length_x_bar, length_y_bar
        complex(dp), intent(out) :: load_bar(:)
        integer, intent(out) :: status

        type(csc_z_t) :: matrix_csc
        type(fortsparse_status_t) :: sparse_status
        type(sparse_direct_factor_t) :: adjoint_factor
        complex(dp), allocatable :: boundary_form(:, :), matrix(:, :)
        complex(dp), allocatable :: matrix_values_bar(:)

        volume_matrix_bar = cmplx(0.0_dp, 0.0_dp, dp)
        sampling_bar = cmplx(0.0_dp, 0.0_dp, dp)
        load_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wave_number_bar = 0.0_dp
        length_x_bar = 0.0_dp
        length_y_bar = 0.0_dp
        call validate_inputs( &
            volume_matrix, sampling, nx, ny, load, state, status)
        if (status /= 0) return
        if (size(state_bar) /= size(state) .or. &
            any(shape(volume_matrix_bar) /= shape(volume_matrix)) .or. &
            any(shape(sampling_bar) /= shape(sampling)) .or. &
            size(load_bar) /= size(load)) then
            status = 1
            return
        end if
        allocate(boundary_form(size(load), size(load)))
        call pullback_planar_maxwell_dtn_form( &
            sampling, nx, ny, wave_number, length_x, length_y, boundary_form, &
            status)
        if (status /= 0) return
        matrix = volume_matrix + boundary_form
        call dense_to_csc(matrix, matrix_csc, sparse_status)
        if (sparse_status%code /= 0) return
        allocate(matrix_values_bar(matrix_csc%nnz))
        call sparse_direct_factor_adjoint_csc( &
            adjoint_factor, size(load), matrix_csc%col_ptr, &
            matrix_csc%row_idx, matrix_csc%val, status)
        if (status /= 0) return
        call sparse_direct_solve_factored_vjp( &
            adjoint_factor, size(load), matrix_csc%col_ptr, &
            matrix_csc%row_idx, state, state_bar, load_bar, matrix_values_bar, &
            status)
        call sparse_direct_free(adjoint_factor)
        if (status /= 0) return
        call csc_values_to_dense( &
            matrix_csc, matrix_values_bar, volume_matrix_bar)
        call pullback_planar_maxwell_dtn_form_vjp( &
            sampling, nx, ny, wave_number, length_x, length_y, &
            volume_matrix_bar, sampling_bar, wave_number_bar, length_x_bar, &
            length_y_bar, status)
    end subroutine solve_planar_maxwell_dtn_system_vjp

    subroutine validate_inputs( &
            volume_matrix, sampling, nx, ny, load, state, status)
        complex(dp), intent(in) :: volume_matrix(:, :), sampling(:, :)
        integer, intent(in) :: nx, ny
        complex(dp), intent(in) :: load(:), state(:)
        integer, intent(out) :: status

        status = 1
        if (size(volume_matrix, 1) /= size(volume_matrix, 2)) return
        if (size(volume_matrix, 1) /= size(load)) return
        if (size(state) /= size(load)) return
        if (size(sampling, 1) /= 2*nx*ny) return
        if (size(sampling, 2) /= size(load)) return
        status = 0
    end subroutine validate_inputs

    subroutine dense_to_csc(matrix, sparse_matrix, sparse_status)
        complex(dp), intent(in) :: matrix(:, :)
        type(csc_z_t), intent(out) :: sparse_matrix
        type(fortsparse_status_t), intent(out) :: sparse_status

        complex(dp), allocatable :: values(:)
        integer, allocatable :: columns(:), rows(:)
        integer :: column, entry, row

        allocate(rows(size(matrix)), columns(size(matrix)), values(size(matrix)))
        entry = 0
        do column = 1, size(matrix, 2)
            do row = 1, size(matrix, 1)
                entry = entry + 1
                rows(entry) = row
                columns(entry) = column
                values(entry) = matrix(row, column)
            end do
        end do
        call csc_from_triplet( &
            size(matrix, 1), size(matrix, 2), rows, columns, values, &
            sparse_matrix, sparse_status)
    end subroutine dense_to_csc

    subroutine csc_values_to_dense(sparse_matrix, values, matrix)
        type(csc_z_t), intent(in) :: sparse_matrix
        complex(dp), intent(in) :: values(:)
        complex(dp), intent(out) :: matrix(:, :)

        integer :: column, entry

        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        do column = 1, sparse_matrix%ncol
            do entry = sparse_matrix%col_ptr(column), &
                    sparse_matrix%col_ptr(column + 1) - 1
                matrix(sparse_matrix%row_idx(entry), column) = values(entry)
            end do
        end do
    end subroutine csc_values_to_dense

end module fortfem_planar_maxwell_dtn_system
