module fortfem_sparse_matrix
    use fortfem_kinds, only: dp
    implicit none
    private

    type :: sparse_matrix_t
        integer :: nrows = 0
        integer :: ncols = 0
        integer, allocatable :: row_ptr(:)
        integer, allocatable :: col_ind(:)
        real(dp), allocatable :: values(:)
    end type sparse_matrix_t

    public :: sparse_matrix_t
    public :: sparse_from_dense
    public :: sparse_to_csc
    public :: spmv

contains

    subroutine sparse_from_dense(A_dense, A_sparse)
        real(dp), intent(in) :: A_dense(:, :)
        type(sparse_matrix_t), intent(out) :: A_sparse

        integer :: nrows, ncols
        integer :: i, j, nnz, idx
        real(dp) :: tol

        nrows = size(A_dense, 1)
        ncols = size(A_dense, 2)

        tol = 1.0e-14_dp

        nnz = 0
        do i = 1, nrows
            do j = 1, ncols
                if (abs(A_dense(i, j)) > tol) then
                    nnz = nnz + 1
                end if
            end do
        end do

        A_sparse%nrows = nrows
        A_sparse%ncols = ncols

        allocate(A_sparse%row_ptr(nrows + 1))
        allocate(A_sparse%col_ind(nnz))
        allocate(A_sparse%values(nnz))

        idx = 1
        A_sparse%row_ptr(1) = 1

        do i = 1, nrows
            do j = 1, ncols
                if (abs(A_dense(i, j)) > tol) then
                    A_sparse%col_ind(idx) = j
                    A_sparse%values(idx) = A_dense(i, j)
                    idx = idx + 1
                end if
            end do
            A_sparse%row_ptr(i + 1) = idx
        end do
    end subroutine sparse_from_dense

    subroutine sparse_to_csc(A_sparse, col_ptr, row_ind, values_csc)
        type(sparse_matrix_t), intent(in) :: A_sparse
        integer, allocatable, intent(out) :: col_ptr(:)
        integer, allocatable, intent(out) :: row_ind(:)
        real(dp), allocatable, intent(out) :: values_csc(:)

        integer :: nrows, ncols, nnz
        integer :: i, j, k, row_start, row_end, idx
        integer, allocatable :: col_counts(:), next(:)

        nrows = A_sparse%nrows
        ncols = A_sparse%ncols
        nnz = size(A_sparse%values)

        allocate (col_ptr(ncols + 1))
        allocate (row_ind(nnz))
        allocate (values_csc(nnz))
        allocate (col_counts(ncols))
        allocate (next(ncols))

        col_counts = 0
        do i = 1, nrows
            row_start = A_sparse%row_ptr(i)
            row_end = A_sparse%row_ptr(i + 1) - 1
            do k = row_start, row_end
                j = A_sparse%col_ind(k)
                col_counts(j) = col_counts(j) + 1
            end do
        end do

        col_ptr(1) = 1
        do j = 1, ncols
            col_ptr(j + 1) = col_ptr(j) + col_counts(j)
        end do

        next(1:ncols) = col_ptr(1:ncols)

        do i = 1, nrows
            row_start = A_sparse%row_ptr(i)
            row_end = A_sparse%row_ptr(i + 1) - 1
            do k = row_start, row_end
                j = A_sparse%col_ind(k)
                idx = next(j)
                row_ind(idx) = i
                values_csc(idx) = A_sparse%values(k)
                next(j) = next(j) + 1
            end do
        end do

        deallocate (col_counts, next)
    end subroutine sparse_to_csc

    subroutine spmv(A, x, y)
        type(sparse_matrix_t), intent(in) :: A
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)

        integer :: i, k, row_start, row_end

        if (size(x) /= A%ncols) then
            error stop "spmv: x has incompatible size"
        end if

        if (size(y) /= A%nrows) then
            error stop "spmv: y has incompatible size"
        end if

        y = 0.0_dp

        do i = 1, A%nrows
            row_start = A%row_ptr(i)
            row_end = A%row_ptr(i + 1) - 1
            do k = row_start, row_end
                y(i) = y(i) + A%values(k) * x(A%col_ind(k))
            end do
        end do
    end subroutine spmv

end module fortfem_sparse_matrix
