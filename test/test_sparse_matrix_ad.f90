program test_sparse_matrix_ad
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_sparse_matrix, only: &
        sparse_from_dense, sparse_matrix_t, spmv, spmv_jvp, spmv_vjp
    implicit none

    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: matrix(3, 3), matrix_dot(3, 3)
    real(dp) :: vector(3), vector_dot(3), product(3), product_dot(3)
    real(dp) :: product_plus(3), product_minus(3), product_dot_fd(3)
    real(dp) :: product_bar(3), vector_bar(3), matrix_values_bar(9)
    real(dp) :: matrix_values_dot(9), matrix_values_plus(9)
    real(dp) :: matrix_values_minus(9), lhs, rhs, relative_error
    type(sparse_matrix_t) :: sparse_matrix, sparse_matrix_plus, sparse_matrix_minus
    integer :: column, entry, row

    matrix = reshape([ &
        4.0_dp, -1.0_dp, 0.5_dp, &
        0.7_dp, 3.0_dp, -0.4_dp, &
        -0.2_dp, 0.8_dp, 2.5_dp], shape(matrix))
    matrix_dot = reshape([ &
        0.12_dp, -0.03_dp, 0.05_dp, &
        -0.07_dp, 0.09_dp, 0.02_dp, &
        0.04_dp, -0.06_dp, 0.08_dp], shape(matrix_dot))
    vector = [0.6_dp, -0.8_dp, 0.35_dp]
    vector_dot = [-0.11_dp, 0.07_dp, 0.09_dp]
    product_bar = [0.31_dp, -0.27_dp, 0.44_dp]

    call sparse_from_dense(matrix, sparse_matrix)
    matrix_values_dot = 0.0_dp
    do column = 1, sparse_matrix%ncol
        do entry = sparse_matrix%col_ptr(column), &
                sparse_matrix%col_ptr(column + 1) - 1
            row = sparse_matrix%row_idx(entry)
            matrix_values_dot(entry) = matrix_dot(row, column)
        end do
    end do
    call spmv(sparse_matrix, vector, product)
    call spmv_jvp( &
        sparse_matrix, vector, matrix_values_dot, vector_dot, product_dot)

    call sparse_from_dense( &
        matrix + step*matrix_dot, sparse_matrix_plus)
    call sparse_from_dense( &
        matrix - step*matrix_dot, sparse_matrix_minus)
    call spmv(sparse_matrix_plus, vector + step*vector_dot, product_plus)
    call spmv(sparse_matrix_minus, vector - step*vector_dot, product_minus)
    product_dot_fd = (product_plus - product_minus)/(2.0_dp*step)
    relative_error = maxval(abs(product_dot - product_dot_fd))/ &
        max(1.0_dp, maxval(abs(product_dot)))
    call check_condition(relative_error < 2.0e-8_dp, &
        "Sparse matrix JVP matches an independent re-evaluation")

    call spmv_vjp( &
        sparse_matrix, vector, product_bar, matrix_values_bar, vector_bar)
    lhs = dot_product(product_bar, product_dot)
    rhs = dot_product(vector_bar, vector_dot) + &
        dot_product(matrix_values_bar, matrix_values_dot)
    call check_condition(abs(lhs - rhs) < &
        2.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Sparse matrix VJP satisfies the adjoint identity")

    matrix_values_plus = sparse_matrix_plus%val
    matrix_values_minus = sparse_matrix_minus%val
    call check_condition(size(matrix_values_plus) == sparse_matrix%nnz .and. &
        size(matrix_values_minus) == sparse_matrix%nnz, &
        "Sparse perturbations preserve the inactive CSC structure")
    call check_summary("Differentiable sparse matrix-vector product")
end program test_sparse_matrix_ad
