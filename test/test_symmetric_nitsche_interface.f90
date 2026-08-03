program test_symmetric_nitsche_interface
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_symmetric_nitsche_interface
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: plus_trace(1, 1) = reshape([1.0_dp], [1, 1])
    real(dp), parameter :: minus_trace(1, 1) = reshape([2.0_dp], [1, 1])
    real(dp), parameter :: plus_flux(1, 1) = reshape([3.0_dp], [1, 1])
    real(dp), parameter :: minus_flux(1, 1) = reshape([4.0_dp], [1, 1])
    real(dp), parameter :: surface_weights(1) = [2.0_dp]
    real(dp), parameter :: penalty = 10.0_dp
    real(dp), parameter :: expected_matrix(2, 2) = reshape([ &
        14.0_dp, -38.0_dp, -38.0_dp, 96.0_dp], [2, 2])
    real(dp) :: matrix(2, 2)
    real(dp), parameter :: bad_weights(2) = [1.0_dp, 2.0_dp]
    type(fortsparse_status_t) :: status

    call assemble_symmetric_nitsche_interface( &
        plus_trace, minus_trace, plus_flux, minus_flux, surface_weights, &
        penalty, matrix, status)
    call check_condition(status%code == 0, &
        "symmetric Nitsche interface accepts compatible traces and fluxes")
    call check_condition(maxval(abs(matrix - expected_matrix)) < 1.0e-14_dp, &
        "symmetric Nitsche block matches the independent consistency oracle")
    call check_condition(maxval(abs(matrix - transpose(matrix))) < 1.0e-14_dp, &
        "symmetric Nitsche interface block is symmetric")

    call assemble_symmetric_nitsche_interface( &
        plus_trace, minus_trace, plus_flux, minus_flux, bad_weights, penalty, &
        matrix, status)
    call check_condition(status%code /= 0, &
        "symmetric Nitsche interface rejects incompatible weights")
    call check_summary("symmetric Nitsche interface")
end program test_symmetric_nitsche_interface
