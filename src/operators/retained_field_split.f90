module fortfem_retained_field_split
    !! Retained block-diagonal field-split solves and derivatives.
    !!
    !! Each field owns one square CSC block and one retained direct factor.
    !! The composition is neutral: coupled off-diagonal blocks remain in a
    !! caller-owned residual or Schur layer, while this module provides the
    !! bounded block-diagonal solve/preconditioner and its fixed-factor JVP/VJP.
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_factor_adjoint_csc, &
        sparse_direct_factor_csc, sparse_direct_factor_t, &
        sparse_direct_factor_transpose_csc, sparse_direct_free, &
        sparse_direct_solve_factored, sparse_direct_solve_factored_jvp, &
        sparse_direct_solve_factored_vjp
    use fortsparse, only: csc_t, csc_z_t, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_NOT_FACTORED, FORTSPARSE_OK
    implicit none
    private

    type, public :: retained_field_split_t
        type(csc_t), allocatable :: matrices(:)
        type(sparse_direct_factor_t), allocatable :: factors(:)
        type(sparse_direct_factor_t), allocatable :: transpose_factors(:)
        integer :: factored_count = 0
    end type retained_field_split_t

    type, public :: retained_complex_field_split_t
        type(csc_z_t), allocatable :: matrices(:)
        type(sparse_direct_factor_t), allocatable :: factors(:)
        type(sparse_direct_factor_t), allocatable :: adjoint_factors(:)
        integer :: factored_count = 0
    end type retained_complex_field_split_t

    public :: factor_retained_field_split
    public :: factor_retained_complex_field_split
    public :: apply_retained_field_split
    public :: apply_retained_complex_field_split
    public :: apply_retained_field_split_jvp
    public :: apply_retained_complex_field_split_jvp
    public :: apply_retained_field_split_vjp
    public :: apply_retained_complex_field_split_vjp
    public :: free_retained_field_split
    public :: free_retained_complex_field_split

contains

    subroutine factor_retained_field_split(matrices, split, status)
        type(csc_t), intent(in) :: matrices(:)
        type(retained_field_split_t), intent(inout) :: split
        integer, intent(out) :: status
        integer :: field

        call free_retained_field_split(split)
        status = FORTSPARSE_INVALID_MATRIX
        if (.not. valid_real_matrices(matrices)) return
        allocate(split%matrices(size(matrices)), split%factors(size(matrices)), &
            split%transpose_factors(size(matrices)))
        do field = 1, size(matrices)
            split%matrices(field) = matrices(field)
            call sparse_direct_factor_csc( &
                split%factors(field), matrices(field)%nrow, matrices(field)%col_ptr, &
                matrices(field)%row_idx, matrices(field)%val, status)
            if (status /= FORTSPARSE_OK) then
                call free_retained_field_split(split)
                return
            end if
            call sparse_direct_factor_transpose_csc( &
                split%transpose_factors(field), matrices(field)%nrow, &
                matrices(field)%col_ptr, matrices(field)%row_idx, matrices(field)%val, status)
            if (status /= FORTSPARSE_OK) then
                call sparse_direct_free(split%factors(field))
                call free_retained_field_split(split)
                return
            end if
            split%factored_count = field
        end do
        call set_status_factored(split%factored_count, status)
    end subroutine factor_retained_field_split

    subroutine factor_retained_complex_field_split(matrices, split, status)
        type(csc_z_t), intent(in) :: matrices(:)
        type(retained_complex_field_split_t), intent(inout) :: split
        integer, intent(out) :: status
        integer :: field

        call free_retained_complex_field_split(split)
        status = FORTSPARSE_INVALID_MATRIX
        if (.not. valid_complex_matrices(matrices)) return
        allocate(split%matrices(size(matrices)), split%factors(size(matrices)), &
            split%adjoint_factors(size(matrices)))
        do field = 1, size(matrices)
            split%matrices(field) = matrices(field)
            call sparse_direct_factor_csc( &
                split%factors(field), matrices(field)%nrow, matrices(field)%col_ptr, &
                matrices(field)%row_idx, matrices(field)%val, status)
            if (status /= FORTSPARSE_OK) then
                call free_retained_complex_field_split(split)
                return
            end if
            call sparse_direct_factor_adjoint_csc( &
                split%adjoint_factors(field), matrices(field)%nrow, &
                matrices(field)%col_ptr, matrices(field)%row_idx, matrices(field)%val, status)
            if (status /= FORTSPARSE_OK) then
                call sparse_direct_free(split%factors(field))
                call free_retained_complex_field_split(split)
                return
            end if
            split%factored_count = field
        end do
        call set_status_factored(split%factored_count, status)
    end subroutine factor_retained_complex_field_split

    subroutine apply_retained_field_split(split, rhs, solution, status)
        type(retained_field_split_t), intent(inout) :: split
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: status
        integer :: field, offset, count

        solution = 0.0_dp
        status = validate_real_split(split, rhs, solution)
        if (status /= FORTSPARSE_OK) return
        offset = 1
        do field = 1, size(split%matrices)
            count = split%matrices(field)%nrow
            call sparse_direct_solve_factored( &
                split%factors(field), rhs(offset:offset + count - 1), &
                solution(offset:offset + count - 1), status)
            if (status /= FORTSPARSE_OK) return
            offset = offset + count
        end do
    end subroutine apply_retained_field_split

    subroutine apply_retained_complex_field_split(split, rhs, solution, status)
        type(retained_complex_field_split_t), intent(inout) :: split
        complex(dp), intent(in) :: rhs(:)
        complex(dp), intent(out) :: solution(:)
        integer, intent(out) :: status
        integer :: field, offset, count

        solution = cmplx(0.0_dp, 0.0_dp, dp)
        status = validate_complex_split(split, rhs, solution)
        if (status /= FORTSPARSE_OK) return
        offset = 1
        do field = 1, size(split%matrices)
            count = split%matrices(field)%nrow
            call sparse_direct_solve_factored( &
                split%factors(field), rhs(offset:offset + count - 1), &
                solution(offset:offset + count - 1), status)
            if (status /= FORTSPARSE_OK) return
            offset = offset + count
        end do
    end subroutine apply_retained_complex_field_split

    subroutine apply_retained_field_split_jvp( &
            split, matrices_dot, solution, rhs_dot, solution_dot, status)
        type(retained_field_split_t), intent(inout) :: split
        type(csc_t), intent(in) :: matrices_dot(:)
        real(dp), intent(in) :: solution(:), rhs_dot(:)
        real(dp), intent(out) :: solution_dot(:)
        integer, intent(out) :: status
        integer :: field, offset, count

        solution_dot = 0.0_dp
        status = validate_real_split_jvp(split, matrices_dot, solution, rhs_dot, solution_dot)
        if (status /= FORTSPARSE_OK) return
        offset = 1
        do field = 1, size(split%matrices)
            count = split%matrices(field)%nrow
            call sparse_direct_solve_factored_jvp( &
                split%factors(field), count, split%matrices(field)%col_ptr, &
                split%matrices(field)%row_idx, matrices_dot(field)%val, &
                solution(offset:offset + count - 1), rhs_dot(offset:offset + count - 1), &
                solution_dot(offset:offset + count - 1), status)
            if (status /= FORTSPARSE_OK) return
            offset = offset + count
        end do
    end subroutine apply_retained_field_split_jvp

    subroutine apply_retained_complex_field_split_jvp( &
            split, matrices_dot, solution, rhs_dot, solution_dot, status)
        type(retained_complex_field_split_t), intent(inout) :: split
        type(csc_z_t), intent(in) :: matrices_dot(:)
        complex(dp), intent(in) :: solution(:), rhs_dot(:)
        complex(dp), intent(out) :: solution_dot(:)
        integer, intent(out) :: status
        integer :: field, offset, count

        solution_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = validate_complex_split_jvp( &
            split, matrices_dot, solution, rhs_dot, solution_dot)
        if (status /= FORTSPARSE_OK) return
        offset = 1
        do field = 1, size(split%matrices)
            count = split%matrices(field)%nrow
            call sparse_direct_solve_factored_jvp( &
                split%factors(field), count, split%matrices(field)%col_ptr, &
                split%matrices(field)%row_idx, matrices_dot(field)%val, &
                solution(offset:offset + count - 1), rhs_dot(offset:offset + count - 1), &
                solution_dot(offset:offset + count - 1), status)
            if (status /= FORTSPARSE_OK) return
            offset = offset + count
        end do
    end subroutine apply_retained_complex_field_split_jvp

    subroutine apply_retained_field_split_vjp( &
            split, solution, solution_bar, rhs_bar, matrices_bar, status)
        type(retained_field_split_t), intent(inout) :: split
        real(dp), intent(in) :: solution(:), solution_bar(:)
        real(dp), intent(out) :: rhs_bar(:)
        type(csc_t), intent(out) :: matrices_bar(:)
        integer, intent(out) :: status
        integer :: field, offset, count

        rhs_bar = 0.0_dp
        status = validate_real_split_vjp(split, solution, solution_bar, rhs_bar, matrices_bar)
        if (status /= FORTSPARSE_OK) return
        offset = 1
        do field = 1, size(split%matrices)
            count = split%matrices(field)%nrow
            call initialize_real_matrix_bar(matrices_bar(field), split%matrices(field))
            call sparse_direct_solve_factored_vjp( &
                split%transpose_factors(field), count, split%matrices(field)%col_ptr, &
                split%matrices(field)%row_idx, solution(offset:offset + count - 1), &
                solution_bar(offset:offset + count - 1), rhs_bar(offset:offset + count - 1), &
                matrices_bar(field)%val, status)
            if (status /= FORTSPARSE_OK) return
            offset = offset + count
        end do
    end subroutine apply_retained_field_split_vjp

    subroutine apply_retained_complex_field_split_vjp( &
            split, solution, solution_bar, rhs_bar, matrices_bar, status)
        type(retained_complex_field_split_t), intent(inout) :: split
        complex(dp), intent(in) :: solution(:), solution_bar(:)
        complex(dp), intent(out) :: rhs_bar(:)
        type(csc_z_t), intent(out) :: matrices_bar(:)
        integer, intent(out) :: status
        integer :: field, offset, count

        rhs_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = validate_complex_split_vjp( &
            split, solution, solution_bar, rhs_bar, matrices_bar)
        if (status /= FORTSPARSE_OK) return
        offset = 1
        do field = 1, size(split%matrices)
            count = split%matrices(field)%nrow
            call initialize_complex_matrix_bar(matrices_bar(field), split%matrices(field))
            call sparse_direct_solve_factored_vjp( &
                split%adjoint_factors(field), count, split%matrices(field)%col_ptr, &
                split%matrices(field)%row_idx, solution(offset:offset + count - 1), &
                solution_bar(offset:offset + count - 1), rhs_bar(offset:offset + count - 1), &
                matrices_bar(field)%val, status)
            if (status /= FORTSPARSE_OK) return
            offset = offset + count
        end do
    end subroutine apply_retained_complex_field_split_vjp

    subroutine free_retained_field_split(split)
        type(retained_field_split_t), intent(inout) :: split
        integer :: field

        if (allocated(split%factors)) then
            do field = 1, min(split%factored_count, size(split%factors))
                call sparse_direct_free(split%factors(field))
                call sparse_direct_free(split%transpose_factors(field))
            end do
        end if
        if (allocated(split%matrices)) deallocate(split%matrices)
        if (allocated(split%factors)) deallocate(split%factors)
        if (allocated(split%transpose_factors)) deallocate(split%transpose_factors)
        split%factored_count = 0
    end subroutine free_retained_field_split

    subroutine free_retained_complex_field_split(split)
        type(retained_complex_field_split_t), intent(inout) :: split
        integer :: field

        if (allocated(split%factors)) then
            do field = 1, min(split%factored_count, size(split%factors))
                call sparse_direct_free(split%factors(field))
                call sparse_direct_free(split%adjoint_factors(field))
            end do
        end if
        if (allocated(split%matrices)) deallocate(split%matrices)
        if (allocated(split%factors)) deallocate(split%factors)
        if (allocated(split%adjoint_factors)) deallocate(split%adjoint_factors)
        split%factored_count = 0
    end subroutine free_retained_complex_field_split

    integer function validate_real_split(split, rhs, solution) result(status)
        type(retained_field_split_t), intent(in) :: split
        real(dp), intent(in) :: rhs(:), solution(:)
        integer :: total

        status = FORTSPARSE_NOT_FACTORED
        if (.not. allocated(split%matrices) .or. .not. allocated(split%factors) .or. &
            .not. allocated(split%transpose_factors) .or. &
            split%factored_count /= size(split%matrices) .or. &
            split%factored_count /= size(split%factors)) return
        total = total_size_real(split%matrices)
        if (total < 1 .or. size(rhs) /= total .or. size(solution) /= total) then
            status = FORTSPARSE_INVALID_MATRIX
            return
        end if
        status = FORTSPARSE_OK
    end function validate_real_split

    integer function validate_complex_split(split, rhs, solution) result(status)
        type(retained_complex_field_split_t), intent(in) :: split
        complex(dp), intent(in) :: rhs(:), solution(:)
        integer :: total

        status = FORTSPARSE_NOT_FACTORED
        if (.not. allocated(split%matrices) .or. .not. allocated(split%factors) .or. &
            .not. allocated(split%adjoint_factors) .or. &
            split%factored_count /= size(split%matrices) .or. &
            split%factored_count /= size(split%factors)) return
        total = total_size_complex(split%matrices)
        if (total < 1 .or. size(rhs) /= total .or. size(solution) /= total) then
            status = FORTSPARSE_INVALID_MATRIX
            return
        end if
        status = FORTSPARSE_OK
    end function validate_complex_split

    integer function validate_real_split_jvp( &
            split, matrices_dot, solution, rhs_dot, solution_dot) result(status)
        type(retained_field_split_t), intent(in) :: split
        type(csc_t), intent(in) :: matrices_dot(:)
        real(dp), intent(in) :: solution(:), rhs_dot(:), solution_dot(:)
        status = validate_real_split(split, rhs_dot, solution_dot)
        if (status /= FORTSPARSE_OK) return
        if (.not. valid_real_matrices(matrices_dot)) then
            status = FORTSPARSE_INVALID_MATRIX
            return
        end if
        if (.not. same_real_patterns(split%matrices, matrices_dot)) status = FORTSPARSE_INVALID_MATRIX
        if (size(solution) /= size(rhs_dot)) status = FORTSPARSE_INVALID_MATRIX
    end function validate_real_split_jvp

    integer function validate_complex_split_jvp( &
            split, matrices_dot, solution, rhs_dot, solution_dot) result(status)
        type(retained_complex_field_split_t), intent(in) :: split
        type(csc_z_t), intent(in) :: matrices_dot(:)
        complex(dp), intent(in) :: solution(:), rhs_dot(:), solution_dot(:)
        status = validate_complex_split(split, rhs_dot, solution_dot)
        if (status /= FORTSPARSE_OK) return
        if (.not. valid_complex_matrices(matrices_dot)) then
            status = FORTSPARSE_INVALID_MATRIX
            return
        end if
        if (.not. same_complex_patterns(split%matrices, matrices_dot)) status = FORTSPARSE_INVALID_MATRIX
        if (size(solution) /= size(rhs_dot)) status = FORTSPARSE_INVALID_MATRIX
    end function validate_complex_split_jvp

    integer function validate_real_split_vjp( &
            split, solution, solution_bar, rhs_bar, matrices_bar) result(status)
        type(retained_field_split_t), intent(in) :: split
        real(dp), intent(in) :: solution(:), solution_bar(:)
        real(dp), intent(in) :: rhs_bar(:)
        type(csc_t), intent(out) :: matrices_bar(:)
        status = validate_real_split(split, solution, rhs_bar)
        if (status /= FORTSPARSE_OK) return
        if (size(solution_bar) /= size(solution) .or. &
            size(matrices_bar) /= size(split%matrices)) status = FORTSPARSE_INVALID_MATRIX
    end function validate_real_split_vjp

    integer function validate_complex_split_vjp( &
            split, solution, solution_bar, rhs_bar, matrices_bar) result(status)
        type(retained_complex_field_split_t), intent(in) :: split
        complex(dp), intent(in) :: solution(:), solution_bar(:)
        complex(dp), intent(in) :: rhs_bar(:)
        type(csc_z_t), intent(out) :: matrices_bar(:)
        status = validate_complex_split(split, solution, rhs_bar)
        if (status /= FORTSPARSE_OK) return
        if (size(solution_bar) /= size(solution) .or. &
            size(matrices_bar) /= size(split%matrices)) status = FORTSPARSE_INVALID_MATRIX
    end function validate_complex_split_vjp

    logical function valid_real_matrices(matrices) result(valid)
        type(csc_t), intent(in) :: matrices(:)
        integer :: field

        valid = size(matrices) > 0
        if (.not. valid) return
        do field = 1, size(matrices)
            if (matrices(field)%nrow < 1 .or. matrices(field)%nrow /= matrices(field)%ncol) then
                valid = .false.
                return
            end if
            if (.not. allocated(matrices(field)%col_ptr) .or. &
                .not. allocated(matrices(field)%row_idx) .or. &
                .not. allocated(matrices(field)%val)) then
                valid = .false.
                return
            end if
            if (size(matrices(field)%col_ptr) /= matrices(field)%ncol + 1 .or. &
                size(matrices(field)%row_idx) /= matrices(field)%nnz .or. &
                size(matrices(field)%val) /= matrices(field)%nnz) then
                valid = .false.
                return
            end if
        end do
    end function valid_real_matrices

    logical function valid_complex_matrices(matrices) result(valid)
        type(csc_z_t), intent(in) :: matrices(:)
        integer :: field

        valid = size(matrices) > 0
        if (.not. valid) return
        do field = 1, size(matrices)
            if (matrices(field)%nrow < 1 .or. matrices(field)%nrow /= matrices(field)%ncol) then
                valid = .false.
                return
            end if
            if (.not. allocated(matrices(field)%col_ptr) .or. &
                .not. allocated(matrices(field)%row_idx) .or. &
                .not. allocated(matrices(field)%val)) then
                valid = .false.
                return
            end if
            if (size(matrices(field)%col_ptr) /= matrices(field)%ncol + 1 .or. &
                size(matrices(field)%row_idx) /= matrices(field)%nnz .or. &
                size(matrices(field)%val) /= matrices(field)%nnz) then
                valid = .false.
                return
            end if
        end do
    end function valid_complex_matrices

    integer function total_size_real(matrices) result(total)
        type(csc_t), intent(in) :: matrices(:)
        integer :: field
        total = 0
        do field = 1, size(matrices)
            total = total + matrices(field)%nrow
        end do
    end function total_size_real

    integer function total_size_complex(matrices) result(total)
        type(csc_z_t), intent(in) :: matrices(:)
        integer :: field
        total = 0
        do field = 1, size(matrices)
            total = total + matrices(field)%nrow
        end do
    end function total_size_complex

    logical function same_real_patterns(left, right) result(same)
        type(csc_t), intent(in) :: left(:), right(:)
        integer :: field
        same = size(left) == size(right)
        if (.not. same) return
        do field = 1, size(left)
            same = same .and. left(field)%nrow == right(field)%nrow .and. &
                left(field)%ncol == right(field)%ncol .and. left(field)%nnz == right(field)%nnz .and. &
                all(left(field)%col_ptr == right(field)%col_ptr) .and. &
                all(left(field)%row_idx == right(field)%row_idx)
        end do
    end function same_real_patterns

    logical function same_complex_patterns(left, right) result(same)
        type(csc_z_t), intent(in) :: left(:), right(:)
        integer :: field
        same = size(left) == size(right)
        if (.not. same) return
        do field = 1, size(left)
            same = same .and. left(field)%nrow == right(field)%nrow .and. &
                left(field)%ncol == right(field)%ncol .and. left(field)%nnz == right(field)%nnz .and. &
                all(left(field)%col_ptr == right(field)%col_ptr) .and. &
                all(left(field)%row_idx == right(field)%row_idx)
        end do
    end function same_complex_patterns

    subroutine initialize_real_matrix_bar(matrix_bar, matrix)
        type(csc_t), intent(out) :: matrix_bar
        type(csc_t), intent(in) :: matrix
        matrix_bar%nrow = matrix%nrow
        matrix_bar%ncol = matrix%ncol
        matrix_bar%nnz = matrix%nnz
        allocate(matrix_bar%col_ptr(size(matrix%col_ptr)), &
            matrix_bar%row_idx(size(matrix%row_idx)), matrix_bar%val(size(matrix%val)))
        matrix_bar%col_ptr = matrix%col_ptr
        matrix_bar%row_idx = matrix%row_idx
        matrix_bar%val = 0.0_dp
    end subroutine initialize_real_matrix_bar

    subroutine initialize_complex_matrix_bar(matrix_bar, matrix)
        type(csc_z_t), intent(out) :: matrix_bar
        type(csc_z_t), intent(in) :: matrix
        matrix_bar%nrow = matrix%nrow
        matrix_bar%ncol = matrix%ncol
        matrix_bar%nnz = matrix%nnz
        allocate(matrix_bar%col_ptr(size(matrix%col_ptr)), &
            matrix_bar%row_idx(size(matrix%row_idx)), matrix_bar%val(size(matrix%val)))
        matrix_bar%col_ptr = matrix%col_ptr
        matrix_bar%row_idx = matrix%row_idx
        matrix_bar%val = cmplx(0.0_dp, 0.0_dp, dp)
    end subroutine initialize_complex_matrix_bar

    subroutine set_status_factored(count, status)
        integer, intent(in) :: count
        integer, intent(out) :: status
        if (count > 0) then
            status = FORTSPARSE_OK
        else
            status = FORTSPARSE_INVALID_MATRIX
        end if
    end subroutine set_status_factored

end module fortfem_retained_field_split
