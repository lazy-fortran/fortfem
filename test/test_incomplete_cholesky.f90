program test_incomplete_cholesky
    use check, only: check_condition, check_summary
    use fortfem_incomplete_cholesky, only: apply_incomplete_cholesky, &
        build_incomplete_cholesky, incomplete_cholesky_factor_t
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: matrix(3, 3) = reshape([ &
        4.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, &
        0.0_dp, 1.0_dp, 2.0_dp], [3, 3])
    real(dp), parameter :: right_hand_side(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: bad_symmetry(2, 2) = reshape([ &
        2.0_dp, 0.0_dp, 1.0_dp, 2.0_dp], [2, 2])
    real(dp), parameter :: bad_definiteness(2, 2) = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp], [2, 2])
    type(incomplete_cholesky_factor_t) :: factor
    real(dp), allocatable :: solution(:)
    integer :: status

    call build_incomplete_cholesky(matrix, factor, status)
    call check_condition(status == 0, "IC(0) accepts a symmetric SPD matrix")
    call apply_incomplete_cholesky( &
        factor, right_hand_side, solution, status)
    call check_condition(status == 0 .and. &
        maxval(abs(matmul(matrix, solution) - right_hand_side)) < 1.0e-12_dp, &
        "IC(0) triangular solves satisfy the independent residual oracle")

    call build_incomplete_cholesky(bad_symmetry, factor, status)
    call check_condition(status /= 0, "IC(0) rejects a nonsymmetric matrix")
    call build_incomplete_cholesky(bad_definiteness, factor, status)
    call check_condition(status /= 0, "IC(0) rejects a non-positive pivot")

    call check_summary("incomplete Cholesky")
end program test_incomplete_cholesky
