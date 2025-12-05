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

