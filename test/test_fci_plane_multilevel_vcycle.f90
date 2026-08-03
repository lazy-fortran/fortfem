program test_fci_plane_multilevel_vcycle
    use check, only: check_condition, check_summary
    use fortfem_feec, only: apply_fci_plane_multilevel_vcycle, &
        apply_fci_plane_multilevel_wcycle
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    type(csc_t) :: operators(3), restrictions(2), prolongations(2)
    type(fortsparse_status_t) :: status
    integer, parameter :: level_offsets(4) = [1, 5, 7, 8]
    real(dp), parameter :: diagonals(7) = [ &
        2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 1.5_dp, 1.5_dp, 1.0_dp]
    real(dp), parameter :: residual(4) = [1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp]
    real(dp), parameter :: fine(4, 4) = reshape([ &
        2.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, &
        -1.0_dp, 2.0_dp, -1.0_dp, 0.0_dp, &
        0.0_dp, -1.0_dp, 2.0_dp, -1.0_dp, &
        0.0_dp, 0.0_dp, -1.0_dp, 2.0_dp], [4, 4])
    real(dp), parameter :: middle(2, 2) = reshape([ &
        1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], [2, 2])
    real(dp), parameter :: coarse(1, 1) = reshape([1.0_dp], [1, 1])
    real(dp), parameter :: restriction_fine(2, 4) = reshape([ &
        0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp], [2, 4])
    real(dp), parameter :: prolongation_fine(4, 2) = reshape([ &
        1.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 1.0_dp], [4, 2])
    real(dp), parameter :: restriction_middle(1, 2) = reshape([0.5_dp, 0.5_dp], [1, 2])
    real(dp), parameter :: prolongation_middle(2, 1) = reshape([1.0_dp, 1.0_dp], [2, 1])
    real(dp) :: correction(4), expected(4)
    integer :: bad_offsets(4)

    call make_matrix( &
        [1, 2, 1, 2, 3, 2, 3, 4, 3, 4], &
        [1, 1, 2, 2, 2, 3, 3, 3, 4, 4], &
        [2.0_dp, -1.0_dp, -1.0_dp, 2.0_dp, -1.0_dp, -1.0_dp, &
        2.0_dp, -1.0_dp, -1.0_dp, 2.0_dp], 4, 4, operators(1), status)
    call make_matrix([1, 2, 1, 2], [1, 1, 2, 2], &
        [1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], 2, 2, operators(2), status)
    call make_matrix([1], [1], [1.0_dp], 1, 1, operators(3), status)
    call make_matrix([1, 1, 2, 2], [1, 2, 3, 4], &
        [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp], 2, 4, restrictions(1), status)
    call make_matrix([1, 2, 3, 4], [1, 1, 2, 2], &
        [1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], 4, 2, prolongations(1), status)
    call make_matrix([1, 1], [1, 2], [0.5_dp, 0.5_dp], &
        1, 2, restrictions(2), status)
    call make_matrix([1, 2], [1, 1], [1.0_dp, 1.0_dp], &
        2, 1, prolongations(2), status)

    call apply_fci_plane_multilevel_vcycle( &
        operators, restrictions, prolongations, level_offsets, diagonals, &
        residual, correction, status)
    call check_condition(status%code == 0, &
        "FCI multilevel plane V-cycle accepts three levels")
    expected = multilevel_oracle( &
        fine, middle, coarse, restriction_fine, prolongation_fine, &
        restriction_middle, prolongation_middle, diagonals, residual)
    call check_condition(maxval(abs(correction - expected)) < 1.0e-13_dp, &
        "FCI multilevel plane V-cycle matches the independent recursive oracle")

    call apply_fci_plane_multilevel_wcycle( &
        operators, restrictions, prolongations, level_offsets, diagonals, &
        residual, correction, status)
    expected = multilevel_w_oracle( &
        fine, middle, coarse, restriction_fine, prolongation_fine, &
        restriction_middle, prolongation_middle, diagonals, residual)
    call check_condition(status%code == 0 .and. &
        maxval(abs(correction - expected)) < 1.0e-13_dp, &
        "FCI multilevel plane W-cycle matches the independent recursive oracle")

    bad_offsets = [1, 5, 6, 8]
    call apply_fci_plane_multilevel_vcycle( &
        operators, restrictions, prolongations, bad_offsets, diagonals, &
        residual, correction, status)
    call check_condition(status%code /= 0, &
        "FCI multilevel plane V-cycle rejects inconsistent level offsets")
    call check_summary("FCI multilevel plane V-cycle")

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

    pure function multilevel_oracle( &
            fine_matrix, middle_matrix, coarse_matrix, restriction_one, &
            prolongation_one, restriction_two, prolongation_two, diagonal, rhs) &
            result(solution)
        real(dp), intent(in) :: fine_matrix(4, 4), middle_matrix(2, 2)
        real(dp), intent(in) :: coarse_matrix(1, 1), restriction_one(2, 4)
        real(dp), intent(in) :: prolongation_one(4, 2), restriction_two(1, 2)
        real(dp), intent(in) :: prolongation_two(2, 1), diagonal(7), rhs(4)
        real(dp) :: solution(4), middle_solution(2), coarse_rhs(1)
        real(dp) :: fine_residual(4), middle_residual(2)

        solution = rhs/diagonal(:4)
        fine_residual = rhs - matmul(fine_matrix, solution)
        middle_solution = matmul(restriction_one, fine_residual)/diagonal(5:6)
        middle_residual = matmul(restriction_one, fine_residual) - &
            matmul(middle_matrix, middle_solution)
        coarse_rhs = matmul(restriction_two, middle_residual)/coarse_matrix(1, 1)
        middle_solution = middle_solution + &
            matmul(prolongation_two, coarse_rhs)
        middle_residual = matmul(restriction_one, fine_residual) - &
            matmul(middle_matrix, middle_solution)
        middle_solution = middle_solution + middle_residual/diagonal(5:6)
        solution = solution + matmul(prolongation_one, middle_solution)
        fine_residual = rhs - matmul(fine_matrix, solution)
        solution = solution + fine_residual/diagonal(:4)
    end function multilevel_oracle

    pure function multilevel_w_oracle( &
            fine_matrix, middle_matrix, coarse_matrix, restriction_one, &
            prolongation_one, restriction_two, prolongation_two, diagonal, rhs) &
            result(solution)
        real(dp), intent(in) :: fine_matrix(4, 4), middle_matrix(2, 2)
        real(dp), intent(in) :: coarse_matrix(1, 1), restriction_one(2, 4)
        real(dp), intent(in) :: prolongation_one(4, 2), restriction_two(1, 2)
        real(dp), intent(in) :: prolongation_two(2, 1), diagonal(7), rhs(4)
        real(dp) :: solution(4), middle_solution(2), coarse_rhs(1)
        real(dp) :: fine_residual(4), middle_rhs(2), middle_residual(2)
        integer :: visit

        solution = rhs/diagonal(:4)
        fine_residual = rhs - matmul(fine_matrix, solution)
        do visit = 1, 2
            middle_rhs = matmul(restriction_one, fine_residual)
            middle_solution = middle_rhs/diagonal(5:6)
            middle_residual = middle_rhs - matmul(middle_matrix, middle_solution)
            coarse_rhs = matmul(restriction_two, middle_residual)/coarse_matrix(1, 1)
            middle_solution = middle_solution + &
                matmul(prolongation_two, coarse_rhs)
            middle_residual = middle_rhs - matmul(middle_matrix, middle_solution)
            middle_solution = middle_solution + middle_residual/diagonal(5:6)
            solution = solution + matmul(prolongation_one, middle_solution)
            if (visit < 2) fine_residual = rhs - matmul(fine_matrix, solution)
        end do
        fine_residual = rhs - matmul(fine_matrix, solution)
        solution = solution + fine_residual/diagonal(:4)
    end function multilevel_w_oracle

end program test_fci_plane_multilevel_vcycle
