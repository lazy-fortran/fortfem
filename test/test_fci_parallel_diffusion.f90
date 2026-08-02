program test_fci_parallel_diffusion
    use check, only: check_condition, check_summary
    use fortfem_feec, only: apply_fci_parallel_diffusion, &
        apply_fci_parallel_diffusion_field_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: forward_map(1, 1, 2) = 1.0_dp
    real(dp), parameter :: backward_map(1, 1, 2) = 1.0_dp
    real(dp), parameter :: line_lengths(1, 2) = reshape( &
        [1.0_dp, 2.0_dp], [1, 2])
    real(dp), parameter :: parallel_coefficient(2) = [2.0_dp, 4.0_dp]
    real(dp), parameter :: canonical_volumes(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: staggered_volumes(2) = [5.0_dp, 7.0_dp]
    real(dp), parameter :: field(3) = [0.0_dp, 1.0_dp, 3.0_dp]
    real(dp), parameter :: expected_diffusion(3) = [ &
        10.0_dp, 2.0_dp, -14.0_dp/3.0_dp]
    real(dp) :: diffusion_field(3)
    real(dp), parameter :: direction(3) = [0.3_dp, -0.5_dp, 0.8_dp]
    real(dp), parameter :: diffusion_field_bar(3) = [0.4_dp, -1.2_dp, 0.7_dp]
    real(dp) :: field_bar(3), directional_action(3)
    real(dp) :: weighted_energy, flux_energy, vjp_left, vjp_right
    type(fortsparse_status_t) :: status

    call apply_fci_parallel_diffusion( &
        forward_map, backward_map, line_lengths, parallel_coefficient, &
        canonical_volumes, staggered_volumes, field, diffusion_field, status)
    call check_condition(status%code == 0, &
        "FCI parallel diffusion accepts positive coefficients and volumes")
    call check_condition(maxval(abs(diffusion_field - expected_diffusion)) < &
        1.0e-14_dp, &
        "FCI parallel diffusion matches the explicit support-operator oracle")

    weighted_energy = dot_product(field*canonical_volumes, diffusion_field)
    flux_energy = -sum(staggered_volumes*parallel_coefficient* &
        [1.0_dp, 1.0_dp]**2)
    call check_condition(abs(weighted_energy - flux_energy) < 1.0e-14_dp, &
        "FCI parallel diffusion satisfies the weighted negative-energy identity")

    call apply_fci_parallel_diffusion( &
        forward_map, backward_map, line_lengths, parallel_coefficient, &
        canonical_volumes, staggered_volumes, direction, directional_action, &
        status)
    call apply_fci_parallel_diffusion_field_vjp( &
        forward_map, backward_map, line_lengths, parallel_coefficient, &
        canonical_volumes, staggered_volumes, diffusion_field_bar, field_bar, &
        status)
    vjp_left = dot_product(diffusion_field_bar, directional_action)
    vjp_right = dot_product(field_bar, direction)
    call check_condition(status%code == 0 .and. &
        abs(vjp_left - vjp_right) < 1.0e-14_dp, &
        "FCI parallel diffusion field VJP satisfies the adjoint identity")

    call apply_fci_parallel_diffusion( &
        forward_map, backward_map, line_lengths, [-1.0_dp, 4.0_dp], &
        canonical_volumes, staggered_volumes, field, diffusion_field, status)
    call check_condition(status%code /= 0, &
        "FCI parallel diffusion rejects a non-positive coefficient")
    call check_summary("FCI parallel diffusion")
end program test_fci_parallel_diffusion
