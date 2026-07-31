program test_fci_parallel_jacobi
    use check, only: check_condition, check_summary
    use fortfem_api, only: apply_fci_parallel_jacobi_preconditioner
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: diagonal(4) = [1.0_dp, 2.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: residual(4) = [2.0_dp, -4.0_dp, 12.0_dp, -5.0_dp]
    real(dp), parameter :: expected_correction(4) = [ &
        2.0_dp, -2.0_dp, 3.0_dp, -1.0_dp]
    real(dp) :: correction(4)
    real(dp), parameter :: bad_diagonal(4) = [1.0_dp, 0.0_dp, 4.0_dp, 5.0_dp]
    type(fortsparse_status_t) :: status

    call apply_fci_parallel_jacobi_preconditioner( &
        diagonal, residual, correction, status)
    call check_condition(status%code == 0, &
        "FCI Jacobi preconditioner accepts a positive diagonal")
    call check_condition(maxval(abs(correction - expected_correction)) < &
        1.0e-14_dp, "FCI Jacobi preconditioner matches diagonal division")

    call apply_fci_parallel_jacobi_preconditioner( &
        bad_diagonal, residual, correction, status)
    call check_condition(status%code /= 0, &
        "FCI Jacobi preconditioner rejects a non-positive diagonal")
    call check_summary("FCI parallel Jacobi preconditioner")
end program test_fci_parallel_jacobi
