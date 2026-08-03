program test_fci_additive_field_split
    use check, only: check_condition, check_summary
    use fortfem_feec, only: apply_fci_additive_field_split_preconditioner
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    type(csc_t) :: fine_operators(2), coarse_operators(2)
    type(csc_t) :: restrictions(2), prolongations(2)
    type(fortsparse_status_t) :: status
    integer, parameter :: plane_offsets(3) = [1, 3, 5]
    real(dp), parameter :: parallel_diagonal(4) = [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: plane_diagonal(4) = [1.5_dp, 1.5_dp, 2.0_dp, 2.0_dp]
    real(dp), parameter :: residual(4) = [6.0_dp, -3.0_dp, 8.0_dp, 10.0_dp]
    real(dp), parameter :: parallel_weight = 0.25_dp
    real(dp), parameter :: plane_weight = 0.75_dp
    real(dp), parameter :: fine(2, 2) = reshape([ &
        1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], [2, 2])
    real(dp), parameter :: coarse(1, 1) = reshape([1.0_dp], [1, 1])
    real(dp), parameter :: restriction(1, 2) = reshape([0.5_dp, 0.5_dp], [1, 2])
    real(dp), parameter :: prolongation(2, 1) = reshape([1.0_dp, 1.0_dp], [2, 1])
    real(dp) :: correction(4), expected(4), plane_correction(4)

    call make_matrix([1, 2, 1, 2], [1, 1, 2, 2], &
        [1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], 2, 2, fine_operators(1), status)
    call make_matrix([1, 2, 1, 2], [1, 1, 2, 2], &
        [1.5_dp, -0.5_dp, -0.5_dp, 1.5_dp], 2, 2, fine_operators(2), status)
    call make_matrix([1], [1], [1.0_dp], 1, 1, coarse_operators(1), status)
    call make_matrix([1], [1], [1.0_dp], 1, 1, coarse_operators(2), status)
    call make_matrix([1, 1], [1, 2], [0.5_dp, 0.5_dp], 1, 2, restrictions(1), status)
    call make_matrix([1, 1], [1, 2], [0.5_dp, 0.5_dp], 1, 2, restrictions(2), status)
    call make_matrix([1, 2], [1, 1], [1.0_dp, 1.0_dp], 2, 1, prolongations(1), status)
    call make_matrix([1, 2], [1, 1], [1.0_dp, 1.0_dp], 2, 1, prolongations(2), status)

    call apply_fci_additive_field_split_preconditioner( &
        parallel_diagonal, plane_diagonal, fine_operators, coarse_operators, &
        restrictions, prolongations, plane_offsets, parallel_weight, &
        plane_weight, residual, correction, status)
    call check_condition(status%code == 0, &
        "FCI additive field split accepts positive diagonal blocks")

    plane_correction(:2) = vcycle( &
        fine, coarse, restriction, prolongation, plane_diagonal(:2), residual(:2))
    plane_correction(3:) = vcycle( &
        fine, coarse, restriction, prolongation, plane_diagonal(3:), residual(3:))
    expected = parallel_weight*residual/parallel_diagonal + &
        plane_weight*plane_correction
    call check_condition(maxval(abs(correction - expected)) < 1.0e-13_dp, &
        "FCI additive field split matches an independent weighted oracle")

    call apply_fci_additive_field_split_preconditioner( &
        [2.0_dp, 0.0_dp, 4.0_dp, 5.0_dp], plane_diagonal, fine_operators, &
        coarse_operators, restrictions, prolongations, plane_offsets, &
        parallel_weight, plane_weight, residual, correction, status)
    call check_condition(status%code /= 0, &
        "FCI additive field split rejects a non-positive parallel diagonal")

    call apply_fci_additive_field_split_preconditioner( &
        parallel_diagonal, plane_diagonal, fine_operators, coarse_operators, &
        restrictions, prolongations, plane_offsets, -0.1_dp, plane_weight, &
        residual, correction, status)
    call check_condition(status%code /= 0, &
        "FCI additive field split rejects a negative partition weight")
    call check_summary("FCI additive field-split preconditioner")

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
            diagonal, right_hand_side) result(correction)
        real(dp), intent(in) :: fine_matrix(2, 2), coarse_matrix(1, 1)
        real(dp), intent(in) :: restriction_matrix(1, 2)
        real(dp), intent(in) :: prolongation_matrix(2, 1)
        real(dp), intent(in) :: diagonal(2), right_hand_side(2)
        real(dp) :: correction(2), fine_residual(2), coarse_rhs(1)

        correction = right_hand_side/diagonal
        fine_residual = right_hand_side - matmul(fine_matrix, correction)
        coarse_rhs = matmul(restriction_matrix, fine_residual)/coarse_matrix(1, 1)
        correction = correction + matmul(prolongation_matrix, coarse_rhs)
        fine_residual = right_hand_side - matmul(fine_matrix, correction)
        correction = correction + fine_residual/diagonal
    end function vcycle

end program test_fci_additive_field_split
