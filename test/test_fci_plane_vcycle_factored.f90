program test_fci_plane_vcycle_factored
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        apply_fci_plane_two_level_vcycle_factored, &
        factor_fci_plane_coarse_operator
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_factor_t, sparse_direct_free
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    type(csc_t) :: fine_operator, coarse_operator
    type(csc_t) :: restriction, prolongation
    type(sparse_direct_factor_t) :: coarse_factor
    type(fortsparse_status_t) :: status
    real(dp), parameter :: fine(3, 3) = reshape([ &
        2.0_dp, -1.0_dp, 0.0_dp, -1.0_dp, 2.0_dp, -1.0_dp, &
        0.0_dp, -1.0_dp, 2.0_dp], [3, 3])
    real(dp), parameter :: coarse(2, 2) = reshape([ &
        1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], [2, 2])
    real(dp), parameter :: restriction_dense(2, 3) = reshape([ &
        0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp], [2, 3])
    real(dp), parameter :: prolongation_dense(3, 2) = reshape([ &
        1.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp, 1.0_dp], [3, 2])
    real(dp), parameter :: diagonal(3) = [2.0_dp, 2.0_dp, 2.0_dp]
    real(dp), parameter :: residual_one(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: residual_two(3) = [-0.7_dp, 0.4_dp, 1.8_dp]
    real(dp) :: correction(3), expected(3)
    call make_matrix( &
        [1, 2, 2, 3, 1, 2, 3, 1], [1, 1, 2, 2, 2, 3, 3, 3], &
        [2.0_dp, -1.0_dp, 2.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, &
        2.0_dp, 0.0_dp], 3, 3, fine_operator, status)
    call make_matrix( &
        [1, 2, 1, 2], [1, 1, 2, 2], &
        [1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], 2, 2, coarse_operator, status)
    call make_matrix( &
        [1, 2, 1, 2], [1, 1, 2, 3], &
        [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp], 2, 3, restriction, status)
    call make_matrix( &
        [1, 2, 2, 3], [1, 1, 2, 2], &
        [1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], 3, 2, prolongation, status)

    call factor_fci_plane_coarse_operator(coarse_operator, coarse_factor, status)
    call check_condition(status%code == 0, &
        "FCI plane coarse factorization accepts the coarse operator")

    call apply_fci_plane_two_level_vcycle_factored( &
        fine_operator, coarse_operator, coarse_factor, restriction, prolongation, &
        diagonal, residual_one, correction, status)
    call check_condition(status%code == 0, &
        "FCI factored plane V-cycle accepts the retained coarse factor")
    expected = vcycle(fine, coarse, restriction_dense, prolongation_dense, &
        diagonal, residual_one)
    call check_condition(maxval(abs(correction - expected)) < 1.0e-13_dp, &
        "FCI factored plane V-cycle matches the independent first oracle")

    call apply_fci_plane_two_level_vcycle_factored( &
        fine_operator, coarse_operator, coarse_factor, restriction, prolongation, &
        diagonal, residual_two, correction, status)
    expected = vcycle(fine, coarse, restriction_dense, prolongation_dense, &
        diagonal, residual_two)
    call check_condition(status%code == 0 .and. &
        maxval(abs(correction - expected)) < 1.0e-13_dp, &
        "FCI factored plane V-cycle reuses the factor for a second right-hand side")

    call sparse_direct_free(coarse_factor)
    call check_summary("FCI retained coarse factor V-cycle")

contains

    subroutine make_matrix(rows, columns, values, row_count, column_count, &
            matrix, matrix_status)
        integer, intent(in) :: rows(:), columns(:), row_count, column_count
        real(dp), intent(in) :: values(:)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: matrix_status

        call csc_from_triplet( &
            row_count, column_count, rows, columns, values, matrix, matrix_status)
    end subroutine make_matrix

    pure function vcycle( &
            fine_matrix, coarse_matrix, restriction_matrix, prolongation_matrix, &
            fine_diagonal, right_hand_side) result(correction_result)
        real(dp), intent(in) :: fine_matrix(3, 3), coarse_matrix(2, 2)
        real(dp), intent(in) :: restriction_matrix(2, 3)
        real(dp), intent(in) :: prolongation_matrix(3, 2)
        real(dp), intent(in) :: fine_diagonal(3), right_hand_side(3)
        real(dp) :: correction_result(3), fine_residual(3), coarse_rhs(2)

        correction_result = right_hand_side/fine_diagonal
        fine_residual = right_hand_side - matmul(fine_matrix, correction_result)
        coarse_rhs = solve_two_by_two( &
            coarse_matrix, matmul(restriction_matrix, fine_residual))
        correction_result = correction_result + &
            matmul(prolongation_matrix, coarse_rhs)
        fine_residual = right_hand_side - matmul(fine_matrix, correction_result)
        correction_result = correction_result + fine_residual/fine_diagonal
    end function vcycle

    pure function solve_two_by_two(matrix, rhs) result(solution)
        real(dp), intent(in) :: matrix(2, 2), rhs(2)
        real(dp) :: solution(2), determinant

        determinant = matrix(1, 1)*matrix(2, 2) - matrix(1, 2)*matrix(2, 1)
        solution = [ &
            (matrix(2, 2)*rhs(1) - matrix(1, 2)*rhs(2))/determinant, &
            (-matrix(2, 1)*rhs(1) + matrix(1, 1)*rhs(2))/determinant]
    end function solve_two_by_two

end program test_fci_plane_vcycle_factored
