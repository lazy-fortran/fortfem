program test_interface_jump_penalty
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_interface_jump_penalty
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: plus_trace(1, 2) = reshape([1.0_dp, 2.0_dp], [1, 2])
    real(dp), parameter :: minus_trace(1, 1) = reshape([3.0_dp], [1, 1])
    real(dp), parameter :: surface_weights(1) = [2.0_dp]
    real(dp), parameter :: penalty = 4.0_dp
    real(dp), parameter :: expected_matrix(3, 3) = reshape([ &
        8.0_dp, 16.0_dp, -24.0_dp, &
        16.0_dp, 32.0_dp, -48.0_dp, &
        -24.0_dp, -48.0_dp, 72.0_dp], [3, 3])
    real(dp) :: penalty_matrix(3, 3), vector(3), energy
    real(dp), parameter :: bad_weights(2) = [1.0_dp, 2.0_dp]
    type(fortsparse_status_t) :: status

    call assemble_interface_jump_penalty( &
        plus_trace, minus_trace, surface_weights, penalty, penalty_matrix, status)
    call check_condition(status%code == 0, &
        "interface penalty accepts compatible plus/minus traces")
    call check_condition(maxval(abs(penalty_matrix - expected_matrix)) < &
        1.0e-14_dp, "interface penalty matches the explicit jump outer-product oracle")
    call check_condition(maxval(abs(penalty_matrix - transpose(penalty_matrix))) < &
        1.0e-14_dp, "interface penalty block is symmetric")
    vector = [0.5_dp, -1.0_dp, 2.0_dp]
    energy = dot_product(vector, matmul(penalty_matrix, vector))
    call check_condition(energy >= -1.0e-14_dp, &
        "interface penalty block is positive semidefinite")

    call assemble_interface_jump_penalty( &
        plus_trace, minus_trace, bad_weights, penalty, penalty_matrix, status)
    call check_condition(status%code /= 0, &
        "interface penalty rejects incompatible surface weights")
    call check_summary("interface jump penalty")
end program test_interface_jump_penalty
