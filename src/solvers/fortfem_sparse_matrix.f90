module fortfem_sparse_matrix
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_matvec, &
        sparse_matrix_t => csc_t, fortsparse_status_t
    implicit none
    private

    public :: sparse_matrix_t
    public :: sparse_from_dense
    public :: spmv
    public :: spmv_jvp
    public :: spmv_vjp

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

    pure subroutine spmv_jvp( &
            matrix, vector, matrix_values_dot, vector_dot, product_dot)
        !! Forward product with the CSC value array and input vector active.
        type(sparse_matrix_t), intent(in) :: matrix
        real(dp), intent(in) :: vector(:), matrix_values_dot(:)
        real(dp), intent(in) :: vector_dot(:)
        real(dp), intent(out) :: product_dot(:)

        integer :: column, entry, row

        if (size(vector) /= matrix%ncol) then
            error stop "spmv_jvp: vector has incompatible size"
        end if
        if (size(vector_dot) /= matrix%ncol) then
            error stop "spmv_jvp: vector_dot has incompatible size"
        end if
        if (size(matrix_values_dot) /= matrix%nnz) then
            error stop "spmv_jvp: matrix value tangent has incompatible size"
        end if
        if (size(product_dot) /= matrix%nrow) then
            error stop "spmv_jvp: product_dot has incompatible size"
        end if

        product_dot = 0.0_dp
        do column = 1, matrix%ncol
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                product_dot(row) = product_dot(row) + &
                    matrix%val(entry)*vector_dot(column) + &
                    matrix_values_dot(entry)*vector(column)
            end do
        end do
    end subroutine spmv_jvp

    pure subroutine spmv_vjp( &
            matrix, vector, product_bar, matrix_values_bar, vector_bar)
        !! Real adjoint of CSC value/vector multiplication.
        type(sparse_matrix_t), intent(in) :: matrix
        real(dp), intent(in) :: vector(:), product_bar(:)
        real(dp), intent(out) :: matrix_values_bar(:), vector_bar(:)

        integer :: column, entry, row

        if (size(vector) /= matrix%ncol) then
            error stop "spmv_vjp: vector has incompatible size"
        end if
        if (size(product_bar) /= matrix%nrow) then
            error stop "spmv_vjp: product_bar has incompatible size"
        end if
        if (size(matrix_values_bar) /= matrix%nnz) then
            error stop "spmv_vjp: matrix value cotangent has incompatible size"
        end if
        if (size(vector_bar) /= matrix%ncol) then
            error stop "spmv_vjp: vector_bar has incompatible size"
        end if

        matrix_values_bar = 0.0_dp
        vector_bar = 0.0_dp
        do column = 1, matrix%ncol
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                matrix_values_bar(entry) = product_bar(row)*vector(column)
                vector_bar(column) = vector_bar(column) + &
                    matrix%val(entry)*product_bar(row)
            end do
        end do
    end subroutine spmv_vjp
end module fortfem_sparse_matrix
