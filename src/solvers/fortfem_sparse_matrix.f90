module fortfem_sparse_matrix
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortsparse, only: csc_from_triplet, csc_matvec, &
        sparse_matrix_t => csc_t, fortsparse_status_t
    implicit none
    private

    public :: sparse_matrix_t
    public :: sparse_from_dense
    public :: csr_matrix_t
    public :: spmv
    public :: spmv_jvp
    public :: spmv_vjp

    integer, parameter :: parallel_row_threshold = 2048

    interface spmv
        module procedure spmv_csc, spmv_csr
    end interface spmv

    interface spmv_jvp
        module procedure spmv_jvp_csc, spmv_jvp_csr
    end interface spmv_jvp

    interface spmv_vjp
        module procedure spmv_vjp_csc, spmv_vjp_csr
    end interface spmv_vjp

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

    subroutine spmv_csc(matrix, vector, product)
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
    end subroutine spmv_csc

    subroutine spmv_jvp_csc( &
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
    end subroutine spmv_jvp_csc

    subroutine spmv_vjp_csc( &
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
    end subroutine spmv_vjp_csc
    subroutine spmv_csr(matrix, vector, product)
        !! Row-parallel CSR product. Each row is summed serially in column order, so
        !! the result is bitwise independent of the thread count. In exchange form
        !! product_i = sum_{j/=i} a_ij (x_j - x_i).
        type(csr_matrix_t), intent(in) :: matrix
        real(dp), intent(in) :: vector(:)
        real(dp), intent(out) :: product(:)

        real(dp) :: accumulator
        integer :: entry, row

        if (size(vector) /= matrix%n) then
            error stop "spmv: vector has incompatible size"
        end if
        if (size(product) /= matrix%n) then
            error stop "spmv: product has incompatible size"
        end if
        !$omp parallel do default(none) shared(matrix, vector, product) &
        !$omp private(row, entry, accumulator) schedule(static) &
        !$omp if(matrix%n >= parallel_row_threshold)
        do row = 1, matrix%n
            accumulator = 0.0_dp
            if (matrix%exchange) then
                do entry = matrix%rowptr(row), matrix%rowptr(row + 1) - 1
                    accumulator = accumulator + matrix%val(entry)* &
                        (vector(matrix%colidx(entry)) - vector(row))
                end do
            else
                do entry = matrix%rowptr(row), matrix%rowptr(row + 1) - 1
                    accumulator = accumulator + &
                        matrix%val(entry)*vector(matrix%colidx(entry))
                end do
            end if
            product(row) = accumulator
        end do
        !$omp end parallel do
    end subroutine spmv_csr

    subroutine spmv_jvp_csr( &
            matrix, vector, matrix_values_dot, vector_dot, product_dot)
        !! Forward product with the CSR value array and input vector active.
        type(csr_matrix_t), intent(in) :: matrix
        real(dp), intent(in) :: vector(:), matrix_values_dot(:)
        real(dp), intent(in) :: vector_dot(:)
        real(dp), intent(out) :: product_dot(:)

        real(dp) :: accumulator, difference, difference_dot
        integer :: column, entry, row

        if (size(vector) /= matrix%n .or. size(vector_dot) /= matrix%n .or. &
            size(product_dot) /= matrix%n) then
            error stop "spmv_jvp: vector has incompatible size"
        end if
        if (size(matrix_values_dot) /= size(matrix%val)) then
            error stop "spmv_jvp: matrix value tangent has incompatible size"
        end if
        !$omp parallel do default(none) &
        !$omp shared(matrix, vector, matrix_values_dot, vector_dot, product_dot) &
        !$omp private(row, entry, column, accumulator, difference, difference_dot) &
        !$omp schedule(static) if(matrix%n >= parallel_row_threshold)
        do row = 1, matrix%n
            accumulator = 0.0_dp
            do entry = matrix%rowptr(row), matrix%rowptr(row + 1) - 1
                column = matrix%colidx(entry)
                difference = vector(column)
                difference_dot = vector_dot(column)
                if (matrix%exchange) then
                    difference = difference - vector(row)
                    difference_dot = difference_dot - vector_dot(row)
                end if
                accumulator = accumulator + matrix%val(entry)*difference_dot + &
                    matrix_values_dot(entry)*difference
            end do
            product_dot(row) = accumulator
        end do
        !$omp end parallel do
    end subroutine spmv_jvp_csr

    subroutine spmv_vjp_csr( &
            matrix, vector, product_bar, matrix_values_bar, vector_bar)
        !! Real adjoint of CSR value/vector multiplication. The symmetric pattern
        !! gives the transposed entry of each entry, so both outputs are written
        !! row-wise without conflicts.
        type(csr_matrix_t), intent(in) :: matrix
        real(dp), intent(in) :: vector(:), product_bar(:)
        real(dp), intent(out) :: matrix_values_bar(:), vector_bar(:)

        real(dp) :: accumulator, row_sum
        integer :: column, entry, row

        if (size(vector) /= matrix%n .or. size(product_bar) /= matrix%n .or. &
            size(vector_bar) /= matrix%n) then
            error stop "spmv_vjp: vector has incompatible size"
        end if
        if (size(matrix_values_bar) /= size(matrix%val)) then
            error stop "spmv_vjp: matrix value cotangent has incompatible size"
        end if
        !$omp parallel do default(none) &
        !$omp shared(matrix, vector, product_bar, matrix_values_bar, vector_bar) &
        !$omp private(row, entry, column, accumulator, row_sum) schedule(static) &
        !$omp if(matrix%n >= parallel_row_threshold)
        do row = 1, matrix%n
            accumulator = 0.0_dp
            row_sum = 0.0_dp
            do entry = matrix%rowptr(row), matrix%rowptr(row + 1) - 1
                column = matrix%colidx(entry)
                if (matrix%exchange) then
                    matrix_values_bar(entry) = product_bar(row)* &
                        (vector(column) - vector(row))
                    row_sum = row_sum + matrix%val(entry)
                else
                    matrix_values_bar(entry) = product_bar(row)*vector(column)
                end if
                accumulator = accumulator + &
                    matrix%val(matrix%tpos(entry))*product_bar(column)
            end do
            if (matrix%exchange) accumulator = accumulator - row_sum*product_bar(row)
            vector_bar(row) = accumulator
        end do
        !$omp end parallel do
    end subroutine spmv_vjp_csr
end module fortfem_sparse_matrix
