program test_fci_parallel_diffusion_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: apply_fci_parallel_diffusion, &
        apply_fci_parallel_diffusion_jvp, apply_fci_parallel_diffusion_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: forward_map(2, 2, 2) = reshape([ &
        0.8_dp, 0.2_dp, 0.1_dp, 0.9_dp, &
        0.7_dp, 0.3_dp, 0.2_dp, 0.8_dp], [2, 2, 2])
    real(dp), parameter :: backward_map(2, 2, 2) = reshape([ &
        0.9_dp, 0.1_dp, 0.2_dp, 0.8_dp, &
        0.85_dp, 0.15_dp, 0.25_dp, 0.75_dp], [2, 2, 2])
    real(dp), parameter :: line_lengths(2, 2) = reshape([ &
        1.2_dp, 1.4_dp, 1.1_dp, 1.3_dp], [2, 2])
    real(dp), parameter :: parallel_coefficient(4) = [ &
        1.5_dp, 1.8_dp, 2.1_dp, 2.4_dp]
    real(dp), parameter :: canonical_volumes(6) = [ &
        1.1_dp, 1.3_dp, 1.5_dp, 1.2_dp, 1.4_dp, 1.6_dp]
    real(dp), parameter :: staggered_volumes(4) = [ &
        0.9_dp, 1.0_dp, 1.1_dp, 1.2_dp]
    real(dp), parameter :: field(6) = [ &
        -0.4_dp, 0.7_dp, 1.1_dp, -0.8_dp, 0.3_dp, 1.5_dp]
    real(dp), parameter :: forward_map_dot(2, 2, 2) = reshape([ &
        0.03_dp, -0.02_dp, 0.01_dp, 0.04_dp, &
        -0.05_dp, 0.06_dp, -0.07_dp, 0.08_dp], [2, 2, 2])
    real(dp), parameter :: backward_map_dot(2, 2, 2) = reshape([ &
        -0.04_dp, 0.05_dp, -0.06_dp, 0.07_dp, &
        0.08_dp, -0.09_dp, 0.10_dp, -0.11_dp], [2, 2, 2])
    real(dp), parameter :: line_lengths_dot(2, 2) = reshape([ &
        0.02_dp, -0.03_dp, 0.04_dp, -0.05_dp], [2, 2])
    real(dp), parameter :: parallel_coefficient_dot(4) = [ &
        -0.06_dp, 0.07_dp, -0.08_dp, 0.09_dp]
    real(dp), parameter :: canonical_volumes_dot(6) = [ &
        0.05_dp, -0.04_dp, 0.03_dp, -0.02_dp, 0.01_dp, 0.06_dp]
    real(dp), parameter :: staggered_volumes_dot(4) = [ &
        -0.03_dp, 0.02_dp, 0.04_dp, -0.01_dp]
    real(dp), parameter :: field_dot(6) = [ &
        0.11_dp, -0.12_dp, 0.13_dp, -0.14_dp, 0.15_dp, -0.16_dp]
    real(dp), parameter :: output_bar(6) = [ &
        0.21_dp, -0.31_dp, 0.41_dp, -0.51_dp, 0.61_dp, -0.71_dp]
    real(dp), parameter :: finite_difference_step = 1.0e-6_dp
    real(dp) :: diffusion_field(6), diffusion_plus(6), diffusion_minus(6)
    real(dp) :: diffusion_dot(6)
    real(dp) :: forward_map_bar(2, 2, 2), backward_map_bar(2, 2, 2)
    real(dp) :: line_lengths_bar(2, 2), parallel_coefficient_bar(4)
    real(dp) :: canonical_volumes_bar(6), staggered_volumes_bar(4), field_bar(6)
    real(dp) :: left_pairing, right_pairing
    type(fortsparse_status_t) :: status

    call apply_fci_parallel_diffusion( &
        forward_map, backward_map, line_lengths, parallel_coefficient, &
        canonical_volumes, staggered_volumes, field, diffusion_field, status)
    call check_condition(status%code == 0, &
        "FCI diffusion AD oracle accepts positive inputs")

    call apply_fci_parallel_diffusion_jvp( &
        forward_map, backward_map, line_lengths, parallel_coefficient, &
        canonical_volumes, staggered_volumes, field, forward_map_dot, &
        backward_map_dot, line_lengths_dot, parallel_coefficient_dot, &
        canonical_volumes_dot, staggered_volumes_dot, field_dot, diffusion_dot, &
        status)
    call apply_fci_parallel_diffusion( &
        forward_map + finite_difference_step*forward_map_dot, &
        backward_map + finite_difference_step*backward_map_dot, &
        line_lengths + finite_difference_step*line_lengths_dot, &
        parallel_coefficient + finite_difference_step*parallel_coefficient_dot, &
        canonical_volumes + finite_difference_step*canonical_volumes_dot, &
        staggered_volumes + finite_difference_step*staggered_volumes_dot, &
        field + finite_difference_step*field_dot, diffusion_plus, status)
    call apply_fci_parallel_diffusion( &
        forward_map - finite_difference_step*forward_map_dot, &
        backward_map - finite_difference_step*backward_map_dot, &
        line_lengths - finite_difference_step*line_lengths_dot, &
        parallel_coefficient - finite_difference_step*parallel_coefficient_dot, &
        canonical_volumes - finite_difference_step*canonical_volumes_dot, &
        staggered_volumes - finite_difference_step*staggered_volumes_dot, &
        field - finite_difference_step*field_dot, diffusion_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(diffusion_dot - (diffusion_plus - diffusion_minus)/ &
        (2.0_dp*finite_difference_step))) < 2.0e-9_dp, &
        "FCI diffusion JVP matches the full finite-difference oracle")

    call apply_fci_parallel_diffusion_vjp( &
        forward_map, backward_map, line_lengths, parallel_coefficient, &
        canonical_volumes, staggered_volumes, field, output_bar, &
        forward_map_bar, backward_map_bar, line_lengths_bar, &
        parallel_coefficient_bar, canonical_volumes_bar, staggered_volumes_bar, &
        field_bar, status)
    left_pairing = dot_product(output_bar, diffusion_dot)
    right_pairing = sum(forward_map_bar*forward_map_dot) + &
        sum(backward_map_bar*backward_map_dot) + &
        sum(line_lengths_bar*line_lengths_dot) + &
        dot_product(parallel_coefficient_bar, parallel_coefficient_dot) + &
        dot_product(canonical_volumes_bar, canonical_volumes_dot) + &
        dot_product(staggered_volumes_bar, staggered_volumes_dot) + &
        dot_product(field_bar, field_dot)
    call check_condition(status%code == 0 .and. &
        abs(left_pairing - right_pairing) < 2.0e-13_dp, &
        "FCI diffusion VJP satisfies the full real dot-product identity")

    call check_summary("FCI parallel diffusion AD")
end program test_fci_parallel_diffusion_ad
