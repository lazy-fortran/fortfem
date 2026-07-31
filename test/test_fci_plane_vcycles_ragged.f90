program test_fci_plane_vcycles_ragged
    use check, only: check_condition, check_summary
    use fortfem_api, only: apply_fci_plane_two_level_vcycles_ragged
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    type(csc_t) :: fine_operators(2), coarse_operators(2)
    type(csc_t) :: restrictions(2), prolongations(2)
    type(fortsparse_status_t) :: status
    integer, parameter :: plane_offsets(3) = [1, 4, 6]
    real(dp), parameter :: diagonal(5) = [2.0_dp, 2.0_dp, 2.0_dp, 1.25_dp, 1.25_dp]
    real(dp), parameter :: residual(5) = [1.0_dp, 2.0_dp, 3.0_dp, 0.4_dp, -1.2_dp]
    real(dp), parameter :: fine_three(3, 3) = reshape([ &
        2.0_dp, -1.0_dp, 0.0_dp, -1.0_dp, 2.0_dp, -1.0_dp, &
        0.0_dp, -1.0_dp, 2.0_dp], [3, 3])
    real(dp), parameter :: coarse_two(2, 2) = reshape([ &
        1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], [2, 2])
    real(dp), parameter :: restriction_three(2, 3) = reshape([ &
        0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp], [2, 3])
    real(dp), parameter :: prolongation_three(3, 2) = reshape([ &
        1.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp, 1.0_dp], [3, 2])
    real(dp), parameter :: fine_two(2, 2) = reshape([ &
        1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], [2, 2])
    real(dp), parameter :: coarse_one(1, 1) = reshape([1.0_dp], [1, 1])
    real(dp), parameter :: restriction_two(1, 2) = reshape([0.5_dp, 0.5_dp], [1, 2])
    real(dp), parameter :: prolongation_two(2, 1) = reshape([1.0_dp, 1.0_dp], [2, 1])
    real(dp) :: correction(5), expected(5)
    integer :: bad_offsets(3)

    call make_matrix( &
        [1, 2, 2, 3, 1, 2, 3, 1], [1, 1, 2, 2, 2, 3, 3, 3], &
        [2.0_dp, -1.0_dp, 2.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, &
        2.0_dp, 0.0_dp], 3, 3, fine_operators(1), status)
    call make_matrix( &
        [1, 2, 1, 2], [1, 1, 2, 2], &
        [1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], 2, 2, &
        coarse_operators(1), status)
    call make_matrix( &
        [1, 2, 1, 2], [1, 1, 2, 3], &
        [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp], 2, 3, &
        restrictions(1), status)
    call make_matrix( &
        [1, 2, 2, 3], [1, 1, 2, 2], &
        [1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], 3, 2, &
        prolongations(1), status)

    call make_matrix( &
        [1, 2, 1, 2], [1, 1, 2, 2], &
        [1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], 2, 2, &
        fine_operators(2), status)
    call make_matrix([1], [1], [1.0_dp], 1, 1, coarse_operators(2), status)
    call make_matrix([1, 1], [1, 2], [0.5_dp, 0.5_dp], 1, 2, &
        restrictions(2), status)
    call make_matrix([1, 2], [1, 1], [1.0_dp, 1.0_dp], 2, 1, &
        prolongations(2), status)

    call apply_fci_plane_two_level_vcycles_ragged( &
        fine_operators, coarse_operators, restrictions, prolongations, &
        plane_offsets, diagonal, residual, correction, status)
    call check_condition(status%code == 0, &
        "FCI ragged plane V-cycle accepts variable plane sizes")

    expected(:3) = vcycle_three( &
        fine_three, coarse_two, restriction_three, prolongation_three, &
        diagonal(:3), residual(:3))
    expected(4:) = vcycle_two( &
        fine_two, coarse_one, restriction_two, prolongation_two, &
        diagonal(4:), residual(4:))
    call check_condition(maxval(abs(correction - expected)) < 1.0e-13_dp, &
        "FCI ragged plane V-cycle matches independent variable-size oracles")

    bad_offsets = [1, 3, 6]
    call apply_fci_plane_two_level_vcycles_ragged( &
        fine_operators, coarse_operators, restrictions, prolongations, &
        bad_offsets, diagonal, residual, correction, status)
    call check_condition(status%code /= 0, &
        "FCI ragged plane V-cycle rejects inconsistent offsets")
    call check_summary("FCI ragged plane two-level V-cycle")

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

    pure function vcycle_three( &
            fine, coarse, restriction, prolongation, diagonal, residual) &
            result(correction)
        real(dp), intent(in) :: fine(3, 3), coarse(2, 2)
        real(dp), intent(in) :: restriction(2, 3), prolongation(3, 2)
        real(dp), intent(in) :: diagonal(3), residual(3)
        real(dp) :: correction(3), fine_residual(3), coarse_rhs(2)

        correction = residual/diagonal
        fine_residual = residual - matmul(fine, correction)
        coarse_rhs = solve_two_by_two(coarse, matmul(restriction, fine_residual))
        correction = correction + matmul(prolongation, coarse_rhs)
        fine_residual = residual - matmul(fine, correction)
        correction = correction + fine_residual/diagonal
    end function vcycle_three

    pure function vcycle_two( &
            fine, coarse, restriction, prolongation, diagonal, residual) &
            result(correction)
        real(dp), intent(in) :: fine(2, 2), coarse(1, 1)
        real(dp), intent(in) :: restriction(1, 2), prolongation(2, 1)
        real(dp), intent(in) :: diagonal(2), residual(2)
        real(dp) :: correction(2), fine_residual(2), coarse_rhs(1)

        correction = residual/diagonal
        fine_residual = residual - matmul(fine, correction)
        coarse_rhs = matmul(restriction, fine_residual)/coarse(1, 1)
        correction = correction + matmul(prolongation, coarse_rhs)
        fine_residual = residual - matmul(fine, correction)
        correction = correction + fine_residual/diagonal
    end function vcycle_two

    pure function solve_two_by_two(matrix, rhs) result(solution)
        real(dp), intent(in) :: matrix(2, 2), rhs(2)
        real(dp) :: solution(2), determinant

        determinant = matrix(1, 1)*matrix(2, 2) - &
            matrix(1, 2)*matrix(2, 1)
        solution = [ &
            (matrix(2, 2)*rhs(1) - matrix(1, 2)*rhs(2))/determinant, &
            (-matrix(2, 1)*rhs(1) + matrix(1, 1)*rhs(2))/determinant]
    end function solve_two_by_two

end program test_fci_plane_vcycles_ragged
