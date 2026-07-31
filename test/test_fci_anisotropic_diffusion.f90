program test_fci_anisotropic_diffusion
    use check, only: check_condition, check_summary
    use fortfem_api, only: apply_fci_anisotropic_diffusion
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
    real(dp), parameter :: field(6) = [0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, &
        2.0_dp, 3.0_dp]
    real(dp), parameter :: expected_diffusion(6) = [ &
        2.0_dp, -1.0_dp, -1.0_dp, -6.0_dp, -7.0_dp, -13.0_dp]
    real(dp) :: diffusion_field(6), weighted_energy
    type(csc_t) :: perpendicular_operators(3)
    type(fortsparse_status_t) :: status

    call make_matrix( &
        [1, 2, 1, 2], [1, 1, 2, 2], [-2.0_dp, 1.0_dp, 1.0_dp, -2.0_dp], &
        perpendicular_operators(1), status)
    call make_matrix( &
        [1, 2], [1, 2], [-1.0_dp, -3.0_dp], perpendicular_operators(2), status)
    call make_matrix( &
        [1, 2], [1, 2], [-3.0_dp, -4.0_dp], perpendicular_operators(3), status)
    call apply_fci_anisotropic_diffusion( &
        perpendicular_operators, forward_map, backward_map, line_lengths, &
        parallel_coefficient, canonical_volumes, staggered_volumes, field, &
        diffusion_field, status)
    call check_condition(status%code == 0, &
        "FCI anisotropic action accepts per-plane and parallel operators")
    call check_condition(maxval(abs(diffusion_field - expected_diffusion)) < &
        1.0e-14_dp, "FCI anisotropic action matches the independent split oracle")
    weighted_energy = dot_product(field, diffusion_field)
    call check_condition(weighted_energy < 0.0_dp, &
        "FCI anisotropic action preserves negative energy for dissipative blocks")

    call apply_fci_anisotropic_diffusion( &
        perpendicular_operators(1:2), forward_map, backward_map, line_lengths, &
        parallel_coefficient, canonical_volumes, staggered_volumes, field, &
        diffusion_field, status)
    call check_condition(status%code /= 0, &
        "FCI anisotropic action rejects a missing plane operator")
    call check_summary("FCI anisotropic diffusion")

contains

    subroutine make_matrix(rows, columns, values, matrix, status)
        integer, intent(in) :: rows(:), columns(:)
        real(dp), intent(in) :: values(:)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        call csc_from_triplet(2, 2, rows, columns, values, matrix, status)
    end subroutine make_matrix

end program test_fci_anisotropic_diffusion
