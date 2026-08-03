program test_fci_plane_vcycle
    use check, only: check_condition, check_summary
    use fortfem_feec, only: apply_fci_plane_two_level_vcycle
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_solve_csc
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    type(csc_t) :: fine_operator, coarse_operator, restriction, prolongation
    type(fortsparse_status_t) :: status
    real(dp), parameter :: diagonal(3) = [2.0_dp, 2.0_dp, 2.0_dp]
    real(dp), parameter :: residual(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp) :: correction(3), expected(3), initial_residual_norm
    real(dp) :: fine_dense(3, 3), coarse_dense(2, 2), restriction_dense(2, 3)
    real(dp) :: prolongation_dense(3, 2), coarse_rhs(2), fine_residual(3)
    real(dp) :: corrected_residual_norm
    real(dp) :: coarse_solution(2), coarse_test_rhs(2)
    integer :: direct_status

    call make_matrix( &
        [1, 2, 2, 3, 1, 2, 3, 1], [1, 1, 2, 2, 2, 3, 3, 3], &
        [2.0_dp, -1.0_dp, 2.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, 2.0_dp, 0.0_dp], &
        3, 3, fine_operator, status)
    call make_matrix( &
        [1, 2, 1, 2], [1, 1, 2, 2], [1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], &
        2, 2, coarse_operator, status)
    call make_matrix( &
        [1, 2, 1, 2], [1, 1, 2, 3], [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp], &
        2, 3, restriction, status)
    call make_matrix( &
        [1, 2, 2, 3], [1, 1, 2, 2], [1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], &
        3, 2, prolongation, status)

    coarse_test_rhs = [1.0_dp, 2.0_dp]
    call sparse_direct_solve_csc( &
        2, coarse_operator%col_ptr, coarse_operator%row_idx, &
        coarse_operator%val, coarse_test_rhs, coarse_solution, direct_status)
    call check_condition(direct_status == 0, &
        "FCI plane coarse matrix has a usable direct solve")

    call apply_fci_plane_two_level_vcycle( &
        fine_operator, coarse_operator, restriction, prolongation, diagonal, &
        residual, correction, status)
    call check_condition(status%code == 0, &
        "FCI plane V-cycle accepts compatible sparse levels")

    fine_dense = reshape([ &
        2.0_dp, -1.0_dp, 0.0_dp, -1.0_dp, 2.0_dp, -1.0_dp, &
        0.0_dp, -1.0_dp, 2.0_dp], [3, 3])
    coarse_dense = reshape([1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], [2, 2])
    restriction_dense = reshape([0.5_dp, 0.5_dp, 0.5_dp, &
        0.0_dp, 0.0_dp, 0.5_dp], [2, 3])
    prolongation_dense = reshape([1.0_dp, 0.5_dp, 0.0_dp, &
        0.0_dp, 0.5_dp, 1.0_dp], [3, 2])
    expected = residual/diagonal
    fine_residual = residual - matmul(fine_dense, expected)
    coarse_rhs = matmul(restriction_dense, fine_residual)
    coarse_rhs = solve_two_by_two(coarse_dense, coarse_rhs)
    expected = expected + matmul(prolongation_dense, coarse_rhs)
    fine_residual = residual - matmul(fine_dense, expected)
    expected = expected + fine_residual/diagonal
    call check_condition(maxval(abs(correction - expected)) < 1.0e-13_dp, &
        "FCI plane V-cycle matches the independent dense two-level oracle")
    initial_residual_norm = sqrt(dot_product(residual, residual))
    corrected_residual_norm = sqrt(dot_product(fine_residual, fine_residual))
    call check_condition(corrected_residual_norm < initial_residual_norm, &
        "FCI plane V-cycle reduces the fine residual")

    call apply_fci_plane_two_level_vcycle( &
        fine_operator, coarse_operator, restriction, prolongation, &
        [-1.0_dp, 2.0_dp, 2.0_dp], residual, correction, status)
    call check_condition(status%code /= 0, &
        "FCI plane V-cycle rejects a non-positive smoother diagonal")
    call check_summary("FCI plane two-level V-cycle")

contains

    subroutine make_matrix(rows, columns, values, row_dimension, &
            column_dimension, matrix, status)
        integer, intent(in) :: rows(:), columns(:), row_dimension
        integer, intent(in) :: column_dimension
        real(dp), intent(in) :: values(:)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        call csc_from_triplet( &
            row_dimension, column_dimension, rows, columns, values, matrix, &
            status)
    end subroutine make_matrix

    pure function solve_two_by_two(matrix, rhs) result(solution)
        real(dp), intent(in) :: matrix(2, 2), rhs(2)
        real(dp) :: solution(2), determinant

        determinant = matrix(1, 1)*matrix(2, 2) - &
            matrix(1, 2)*matrix(2, 1)
        solution = [ &
            (matrix(2, 2)*rhs(1) - matrix(1, 2)*rhs(2))/determinant, &
            (-matrix(2, 1)*rhs(1) + matrix(1, 1)*rhs(2))/determinant]
    end function solve_two_by_two

end program test_fci_plane_vcycle
