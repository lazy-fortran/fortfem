program test_fci_anisotropic_diagonal
    use check, only: check_condition, check_summary
    use fortfem_api, only: compute_fci_anisotropic_diffusion_diagonal
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    real(dp), parameter :: forward_map(2, 2, 2) = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2, 2])
    real(dp), parameter :: backward_map(2, 2, 2) = forward_map
    real(dp), parameter :: line_lengths(2, 2) = 1.0_dp
    real(dp), parameter :: parallel_coefficient(4) = 1.0_dp
    real(dp), parameter :: canonical_volumes(6) = 1.0_dp
    real(dp), parameter :: staggered_volumes(4) = 1.0_dp
    real(dp), parameter :: expected_diagonal(6) = [ &
        3.0_dp, 3.0_dp, 3.0_dp, 5.0_dp, 4.0_dp, 5.0_dp]
    real(dp) :: diagonal(6)
    type(csc_t) :: perpendicular_operators(3)
    type(fortsparse_status_t) :: status

    call make_matrix( &
        [1, 2, 1, 2], [1, 1, 2, 2], [-2.0_dp, 1.0_dp, 1.0_dp, -2.0_dp], &
        perpendicular_operators(1), status)
    call make_matrix( &
        [1, 2], [1, 2], [-1.0_dp, -3.0_dp], perpendicular_operators(2), &
        status)
    call make_matrix( &
        [1, 2], [1, 2], [-3.0_dp, -4.0_dp], perpendicular_operators(3), &
        status)

    call compute_fci_anisotropic_diffusion_diagonal( &
        perpendicular_operators, forward_map, backward_map, line_lengths, &
        parallel_coefficient, canonical_volumes, staggered_volumes, diagonal, &
        status)
    call check_condition(status%code == 0, &
        "FCI anisotropic diagonal accepts dissipative plane blocks")
    call check_condition(maxval(abs(diagonal - expected_diagonal)) < &
        1.0e-14_dp, "FCI anisotropic diagonal matches the independent oracle")
    call check_condition(all(diagonal > 0.0_dp), &
        "FCI anisotropic diagonal is positive")

    call make_matrix( &
        [1, 2], [1, 2], [1.0_dp, 1.0_dp], perpendicular_operators(1), status)
    call compute_fci_anisotropic_diffusion_diagonal( &
        perpendicular_operators, forward_map, backward_map, line_lengths, &
        parallel_coefficient, canonical_volumes, staggered_volumes, diagonal, &
        status)
    call check_condition(status%code /= 0, &
        "FCI anisotropic diagonal rejects a non-positive combined diagonal")
    call check_summary("FCI anisotropic diffusion diagonal")

contains

    subroutine make_matrix(rows, columns, values, matrix, status)
        integer, intent(in) :: rows(:), columns(:)
        real(dp), intent(in) :: values(:)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        call csc_from_triplet(2, 2, rows, columns, values, matrix, status)
    end subroutine make_matrix

end program test_fci_anisotropic_diagonal
