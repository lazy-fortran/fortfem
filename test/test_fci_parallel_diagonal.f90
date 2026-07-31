program test_fci_parallel_diagonal
    use check, only: check_condition, check_summary
    use fortfem_api, only: compute_fci_parallel_diffusion_diagonal
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: forward_map(1, 1, 2) = 1.0_dp
    real(dp), parameter :: backward_map(1, 1, 2) = 1.0_dp
    real(dp), parameter :: line_lengths(1, 2) = reshape([ &
        1.0_dp, 2.0_dp], [1, 2])
    real(dp), parameter :: parallel_coefficient(2) = [2.0_dp, 4.0_dp]
    real(dp), parameter :: canonical_volumes(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: staggered_volumes(2) = [5.0_dp, 7.0_dp]
    real(dp), parameter :: expected_diagonal(3) = [ &
        10.0_dp, 8.5_dp, 7.0_dp/3.0_dp]
    real(dp) :: diagonal(3)
    real(dp), parameter :: bad_lengths(1, 2) = reshape([ &
        1.0_dp, -2.0_dp], [1, 2])
    type(fortsparse_status_t) :: status

    call compute_fci_parallel_diffusion_diagonal( &
        forward_map, backward_map, line_lengths, parallel_coefficient, &
        canonical_volumes, staggered_volumes, diagonal, status)
    call check_condition(status%code == 0, &
        "FCI diffusion diagonal accepts positive support data")
    call check_condition(maxval(abs(diagonal - expected_diagonal)) < &
        1.0e-14_dp, "FCI diffusion diagonal matches the explicit Q-squared oracle")
    call check_condition(all(diagonal > 0.0_dp), &
        "FCI diffusion diagonal is positive for nonzero support weights")

    call compute_fci_parallel_diffusion_diagonal( &
        forward_map, backward_map, bad_lengths, parallel_coefficient, &
        canonical_volumes, staggered_volumes, diagonal, status)
    call check_condition(status%code /= 0, &
        "FCI diffusion diagonal rejects non-positive line lengths")
    call check_summary("FCI parallel diffusion diagonal")
end program test_fci_parallel_diagonal
