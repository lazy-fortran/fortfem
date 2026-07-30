module fortfem_sparse_matrix
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_matvec, &
        sparse_matrix_t => csc_t, fortsparse_status_t
    implicit none
    private

    public :: sparse_matrix_t
    public :: sparse_from_dense
    public :: spmv

contains

    subroutine sparse_from_dense(dense_matrix, sparse_matrix)
        real(dp), intent(in) :: dense_matrix(:, :)
        type(sparse_matrix_t), intent(out) :: sparse_matrix

        type(fortsparse_status_t) :: status
        integer, allocatable :: columns(:), rows(:)
        real(dp), allocatable :: values(:)
        integer :: column, entry, nonzero_count, row

        nonzero_count = count(abs(dense_matrix) > 1.0e-14_dp)
        allocate(rows(nonzero_count), columns(nonzero_count))
        allocate(values(nonzero_count))
        entry = 0
        do column = 1, size(dense_matrix, 2)
            do row = 1, size(dense_matrix, 1)
                if (abs(dense_matrix(row, column)) <= 1.0e-14_dp) cycle
                entry = entry + 1
                rows(entry) = row
                columns(entry) = column
                values(entry) = dense_matrix(row, column)
            end do
        end do
        call csc_from_triplet( &
            size(dense_matrix, 1), size(dense_matrix, 2), rows, columns, &
            values, sparse_matrix, status)
        if (status%code /= 0) then
            error stop "sparse_from_dense: FortSparse CSC construction failed"
        end if
    end subroutine sparse_from_dense

    pure subroutine spmv(matrix, vector, product)
        type(sparse_matrix_t), intent(in) :: matrix
        real(dp), intent(in) :: vector(:)
        real(dp), intent(out) :: product(:)

        if (size(vector) /= matrix%ncol) then
            error stop "spmv: vector has incompatible size"
        end if
        if (size(product) /= matrix%nrow) then
            error stop "spmv: product has incompatible size"
        end if
        product = csc_matvec(matrix, vector)
    end subroutine spmv
end module fortfem_sparse_matrix
