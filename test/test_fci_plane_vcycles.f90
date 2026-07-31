program test_fci_plane_vcycles
    use check, only: check_condition, check_summary
    use fortfem_api, only: apply_fci_plane_two_level_vcycles
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    type(csc_t) :: fine_operators(2), coarse_operators(2)
    type(csc_t) :: restrictions(2), prolongations(2)
    type(fortsparse_status_t) :: status
    real(dp), parameter :: fine_dense(3, 3) = reshape([ &
        2.0_dp, -1.0_dp, 0.0_dp, -1.0_dp, 2.0_dp, -1.0_dp, &
        0.0_dp, -1.0_dp, 2.0_dp], [3, 3])
    real(dp), parameter :: coarse_dense(2, 2) = reshape([ &
        1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], [2, 2])
    real(dp), parameter :: restriction_dense(2, 3) = reshape([ &
        0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp], [2, 3])
    real(dp), parameter :: prolongation_dense(3, 2) = reshape([ &
        1.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp, 1.0_dp], [3, 2])
    real(dp), parameter :: diagonal(3, 2) = reshape([ &
        2.0_dp, 2.0_dp, 2.0_dp, 1.5_dp, 2.5_dp, 3.0_dp], [3, 2])
    real(dp), parameter :: residual(6) = [1.0_dp, 2.0_dp, 3.0_dp, -1.0_dp, 0.5_dp, 2.0_dp]
    real(dp) :: correction(6), expected(6), bad_diagonal(2, 2)
    integer :: plane

    do plane = 1, 2
        call make_matrix( &
            [1, 2, 2, 3, 1, 2, 3, 1], [1, 1, 2, 2, 2, 3, 3, 3], &
            [2.0_dp, -1.0_dp, 2.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, &
            2.0_dp, 0.0_dp], 3, 3, fine_operators(plane), status)
        call make_matrix( &
            [1, 2, 1, 2], [1, 1, 2, 2], &
            [1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], 2, 2, &
            coarse_operators(plane), status)
        call make_matrix( &
            [1, 2, 1, 2], [1, 1, 2, 3], &
            [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp], 2, 3, &
            restrictions(plane), status)
        call make_matrix( &
            [1, 2, 2, 3], [1, 1, 2, 2], &
            [1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], 3, 2, &
            prolongations(plane), status)
    end do

    call apply_fci_plane_two_level_vcycles( &
        fine_operators, coarse_operators, restrictions, prolongations, &
        diagonal, residual, correction, status)
    call check_condition(status%code == 0, &
        "FCI batched plane V-cycle accepts compatible levels")

    expected(:3) = vcycle_oracle( &
        fine_dense, coarse_dense, restriction_dense, prolongation_dense, &
        diagonal(:, 1), residual(:3))
    expected(4:) = vcycle_oracle( &
        fine_dense, coarse_dense, restriction_dense, prolongation_dense, &
        diagonal(:, 2), residual(4:))
    call check_condition(maxval(abs(correction - expected)) < 1.0e-13_dp, &
        "FCI batched plane V-cycle matches independent block oracles")

    bad_diagonal = 0.0_dp
    call apply_fci_plane_two_level_vcycles( &
        fine_operators, coarse_operators, restrictions, prolongations, &
        bad_diagonal, residual, correction, status)
    call check_condition(status%code /= 0, &
        "FCI batched plane V-cycle rejects incompatible diagonal storage")
    call check_summary("FCI batched plane two-level V-cycle")

contains

    subroutine make_matrix(rows, columns, values, row_dimension, &
            column_dimension, matrix, matrix_status)
        integer, intent(in) :: rows(:), columns(:), row_dimension
        integer, intent(in) :: column_dimension
        real(dp), intent(in) :: values(:)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: matrix_status

        call csc_from_triplet( &
            row_dimension, column_dimension, rows, columns, values, matrix, &
            matrix_status)
    end subroutine make_matrix

    pure function vcycle_oracle( &
            fine, coarse, restriction, prolongation, diagonal, residual) &
            result(correction)
        real(dp), intent(in) :: fine(3, 3), coarse(2, 2)
        real(dp), intent(in) :: restriction(2, 3), prolongation(3, 2)
        real(dp), intent(in) :: diagonal(3), residual(3)
        real(dp) :: correction(3), fine_residual(3), coarse_rhs(2)

        correction = residual/diagonal
        fine_residual = residual - matmul(fine, correction)
        coarse_rhs = matmul(restriction, fine_residual)
        coarse_rhs = solve_two_by_two(coarse, coarse_rhs)
        correction = correction + matmul(prolongation, coarse_rhs)
        fine_residual = residual - matmul(fine, correction)
        correction = correction + fine_residual/diagonal
    end function vcycle_oracle

    pure function solve_two_by_two(matrix, rhs) result(solution)
        real(dp), intent(in) :: matrix(2, 2), rhs(2)
        real(dp) :: solution(2), determinant

        determinant = matrix(1, 1)*matrix(2, 2) - &
            matrix(1, 2)*matrix(2, 1)
        solution = [ &
            (matrix(2, 2)*rhs(1) - matrix(1, 2)*rhs(2))/determinant, &
            (-matrix(2, 1)*rhs(1) + matrix(1, 1)*rhs(2))/determinant]
    end function solve_two_by_two

end program test_fci_plane_vcycles
