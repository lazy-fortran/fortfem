program test_mortar_trace_coupling
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_mortar_trace_coupling
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: test_trace(2, 2) = reshape([ &
        1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2, 2])
    real(dp), parameter :: trial_trace(2, 1) = reshape([5.0_dp, 6.0_dp], [2, 1])
    real(dp), parameter :: surface_weights(2) = [2.0_dp, 1.0_dp]
    real(dp), parameter :: expected_matrix(2, 1) = reshape([22.0_dp, 54.0_dp], [2, 1])
    real(dp) :: matrix(2, 1)
    real(dp), parameter :: bad_weights(1) = [1.0_dp]
    type(fortsparse_status_t) :: status

    call assemble_mortar_trace_coupling( &
        test_trace, trial_trace, surface_weights, matrix, status)
    call check_condition(status%code == 0, &
        "mortar trace coupling accepts independent test/trial traces")
    call check_condition(maxval(abs(matrix - expected_matrix)) < 1.0e-14_dp, &
        "mortar trace coupling matches the cross-mass oracle")

    call assemble_mortar_trace_coupling( &
        test_trace, trial_trace, bad_weights, matrix, status)
    call check_condition(status%code /= 0, &
        "mortar trace coupling rejects incompatible quadrature sizes")
    call check_summary("mortar trace coupling")
end program test_mortar_trace_coupling
