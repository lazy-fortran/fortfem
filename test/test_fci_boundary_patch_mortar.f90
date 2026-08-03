program test_fci_boundary_patch_mortar
    use check, only: check_condition, check_summary
    use fortfem_fci_boundary_patch_mortar, only: assemble_fci_boundary_patch_mortar, &
        assemble_fci_boundary_patch_mortar_jvp, &
        assemble_fci_boundary_patch_mortar_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: tolerance = 2.0e-13_dp
    real(dp) :: background_trace(3, 2), patch_trace(3, 3)
    real(dp) :: background_trace_reversed(3, 2), patch_trace_reversed(3, 3)
    real(dp) :: rank_deficient_trace(3, 2)
    real(dp) :: overlap_weights(3), reversed_weights(3), zero_weights(3)
    real(dp) :: cross_mass(2, 3), cross_mass_reversed(2, 3)
    real(dp) :: background_mass(2), patch_mass(3)
    real(dp) :: background_from_patch(2, 3), patch_from_background(3, 2)
    real(dp) :: background_mass_reversed(2), patch_mass_reversed(3)
    real(dp) :: background_transfer_reversed(2, 3)
    real(dp) :: patch_transfer_reversed(3, 2)
    real(dp) :: expected_cross_mass(2, 3)
    real(dp) :: expected_background_mass(2), expected_patch_mass(3)
    real(dp) :: background_trace_dot(3, 2), patch_trace_dot(3, 3)
    real(dp) :: overlap_weights_dot(3)
    real(dp) :: cross_mass_dot(2, 3), background_mass_dot(2), patch_mass_dot(3)
    real(dp) :: background_from_patch_dot(2, 3)
    real(dp) :: patch_from_background_dot(3, 2), overlap_measure_dot
    real(dp) :: cross_mass_plus(2, 3), cross_mass_minus(2, 3)
    real(dp) :: background_mass_plus(2), background_mass_minus(2)
    real(dp) :: patch_mass_plus(3), patch_mass_minus(3)
    real(dp) :: background_from_patch_plus(2, 3)
    real(dp) :: background_from_patch_minus(2, 3)
    real(dp) :: patch_from_background_plus(3, 2)
    real(dp) :: patch_from_background_minus(3, 2)
    real(dp) :: overlap_measure_plus, overlap_measure_minus
    real(dp) :: background_trace_bar(3, 2), patch_trace_bar(3, 3)
    real(dp) :: overlap_weights_bar(3)
    real(dp) :: cross_mass_bar(2, 3), background_mass_bar(2), patch_mass_bar(3)
    real(dp) :: background_from_patch_bar(2, 3)
    real(dp) :: patch_from_background_bar(3, 2), overlap_measure_bar
    real(dp) :: finite_difference(2, 3), finite_difference_patch(3, 2)
    real(dp) :: forward_pairing, reverse_pairing
    real(dp) :: overlap_measure, overlap_measure_reversed
    real(dp) :: state_background(2), state_patch(3)
    real(dp) :: left_pairing, right_pairing
    integer :: ownership(3), reversed_ownership(3)
    type(fortsparse_status_t) :: status

    background_trace(:, 1) = [1.0_dp, 0.5_dp, 0.0_dp]
    background_trace(:, 2) = [0.0_dp, 0.5_dp, 1.0_dp]
    patch_trace(:, 1) = [1.0_dp, 0.5_dp, 0.0_dp]
    patch_trace(:, 2) = [0.0_dp, 0.5_dp, 0.5_dp]
    patch_trace(:, 3) = [0.0_dp, 0.0_dp, 0.5_dp]
    overlap_weights = [1.0_dp, 2.0_dp, 1.0_dp]
    ownership = [1, 1, 1]

    expected_cross_mass = 0.0_dp
    expected_cross_mass(1, 1) = 1.5_dp
    expected_cross_mass(1, 2) = 0.5_dp
    expected_cross_mass(1, 3) = 0.0_dp
    expected_cross_mass(2, 1) = 0.5_dp
    expected_cross_mass(2, 2) = 1.0_dp
    expected_cross_mass(2, 3) = 0.5_dp
    expected_background_mass = sum(expected_cross_mass, dim=2)
    expected_patch_mass = sum(expected_cross_mass, dim=1)

    call assemble_fci_boundary_patch_mortar( &
        background_trace, patch_trace, overlap_weights, ownership, cross_mass, &
        background_mass, patch_mass, background_from_patch, &
        patch_from_background, overlap_measure, status)
    call check_condition(status%code == 0, &
        "FCI boundary-patch mortar accepts a nonmatching overlap")
    call check_condition(maxval(abs(cross_mass - expected_cross_mass)) < tolerance, &
        "FCI boundary-patch cross-mass matches an independent oracle")
    call check_condition(maxval(abs(background_mass - expected_background_mass)) < &
        tolerance, "background lumped measure matches the cross-mass row sum")
    call check_condition(maxval(abs(patch_mass - expected_patch_mass)) < tolerance, &
        "patch lumped measure matches the cross-mass column sum")
    call check_condition(abs(overlap_measure - 4.0_dp) < tolerance, &
        "FCI boundary-patch mortar reports the analytical overlap measure")

    call check_condition(maxval(abs(matmul(background_from_patch, &
        [1.0_dp, 1.0_dp, 1.0_dp]) - [1.0_dp, 1.0_dp])) < tolerance, &
        "background transfer reproduces constants")
    call check_condition(maxval(abs(matmul(patch_from_background, &
        [1.0_dp, 1.0_dp]) - [1.0_dp, 1.0_dp, 1.0_dp])) < tolerance, &
        "patch transfer reproduces constants")

    state_background = [0.7_dp, -0.4_dp]
    state_patch = [1.2_dp, -0.3_dp, 0.8_dp]
    left_pairing = dot_product(state_background, background_mass* &
        matmul(background_from_patch, state_patch))
    right_pairing = dot_product(matmul(patch_from_background, state_background), &
        patch_mass*state_patch)
    call check_condition(abs(left_pairing - right_pairing) < tolerance, &
        "paired transfer operators satisfy the weighted adjoint identity")

    background_trace_dot(:, 1) = [0.10_dp, 0.03_dp, -0.04_dp]
    background_trace_dot(:, 2) = -background_trace_dot(:, 1)
    patch_trace_dot(:, 1) = [0.06_dp, -0.02_dp, 0.05_dp]
    patch_trace_dot(:, 2) = [-0.01_dp, 0.04_dp, -0.03_dp]
    patch_trace_dot(:, 3) = -patch_trace_dot(:, 1) - patch_trace_dot(:, 2)
    overlap_weights_dot = [0.20_dp, -0.10_dp, 0.15_dp]
    call assemble_fci_boundary_patch_mortar_jvp( &
        background_trace, patch_trace, overlap_weights, ownership, &
        background_trace_dot, patch_trace_dot, overlap_weights_dot, &
        cross_mass_dot, background_mass_dot, patch_mass_dot, &
        background_from_patch_dot, patch_from_background_dot, &
        overlap_measure_dot, status)
    call check_condition(status%code == 0, &
        "FCI boundary-patch mortar fixed-topology JVP succeeds")

    call assemble_fci_boundary_patch_mortar( &
        background_trace + 1.0e-6_dp*background_trace_dot, &
        patch_trace + 1.0e-6_dp*patch_trace_dot, &
        overlap_weights + 1.0e-6_dp*overlap_weights_dot, ownership, &
        cross_mass_plus, background_mass_plus, patch_mass_plus, &
        background_from_patch_plus, patch_from_background_plus, &
        overlap_measure_plus, status)
    call assemble_fci_boundary_patch_mortar( &
        background_trace - 1.0e-6_dp*background_trace_dot, &
        patch_trace - 1.0e-6_dp*patch_trace_dot, &
        overlap_weights - 1.0e-6_dp*overlap_weights_dot, ownership, &
        cross_mass_minus, background_mass_minus, patch_mass_minus, &
        background_from_patch_minus, patch_from_background_minus, &
        overlap_measure_minus, status)
    finite_difference = (cross_mass_plus - cross_mass_minus)/(2.0e-6_dp)
    call check_condition(maxval(abs(cross_mass_dot - finite_difference)) < 5.0e-10_dp, &
        "FCI boundary-patch cross-mass JVP matches finite differences")
    finite_difference = (background_from_patch_plus - &
        background_from_patch_minus)/(2.0e-6_dp)
    call check_condition(maxval(abs(background_from_patch_dot - finite_difference)) < &
        5.0e-10_dp, "background transfer JVP matches finite differences")
    finite_difference_patch = (patch_from_background_plus - &
        patch_from_background_minus)/(2.0e-6_dp)
    call check_condition(maxval(abs(patch_from_background_dot - &
        finite_difference_patch)) < &
        5.0e-10_dp, "patch transfer JVP matches finite differences")
    call check_condition(abs(overlap_measure_dot - &
        (overlap_measure_plus - overlap_measure_minus)/(2.0e-6_dp)) < 5.0e-10_dp, &
        "overlap-measure JVP matches finite differences")

    cross_mass_bar = reshape([0.3_dp, -0.2_dp, 0.1_dp, -0.4_dp, &
        0.5_dp, -0.6_dp], [2, 3])
    background_mass_bar = [0.7_dp, -0.8_dp]
    patch_mass_bar = [0.2_dp, -0.3_dp, 0.4_dp]
    background_from_patch_bar = reshape([0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp, &
        0.5_dp, -0.6_dp], [2, 3])
    patch_from_background_bar = reshape([-0.2_dp, 0.3_dp, -0.4_dp, 0.5_dp, &
        -0.6_dp, 0.7_dp], [3, 2])
    overlap_measure_bar = 0.9_dp
    call assemble_fci_boundary_patch_mortar_vjp( &
        background_trace, patch_trace, overlap_weights, ownership, &
        cross_mass_bar, background_mass_bar, patch_mass_bar, &
        background_from_patch_bar, patch_from_background_bar, &
        overlap_measure_bar, background_trace_bar, patch_trace_bar, &
        overlap_weights_bar, status)
    call check_condition(status%code == 0, &
        "FCI boundary-patch mortar fixed-topology VJP succeeds")
    forward_pairing = sum(cross_mass_bar*cross_mass_dot) + &
        sum(background_mass_bar*background_mass_dot) + &
        sum(patch_mass_bar*patch_mass_dot) + &
        sum(background_from_patch_bar*background_from_patch_dot) + &
        sum(patch_from_background_bar*patch_from_background_dot) + &
        overlap_measure_bar*overlap_measure_dot
    reverse_pairing = sum(background_trace_bar*background_trace_dot) + &
        sum(patch_trace_bar*patch_trace_dot) + &
        sum(overlap_weights_bar*overlap_weights_dot)
    call check_condition(abs(forward_pairing - reverse_pairing) < 5.0e-11_dp, &
        "FCI boundary-patch mortar VJP satisfies the real dot-product identity")

    background_trace_reversed = background_trace(3:1:-1, :)
    patch_trace_reversed = patch_trace(3:1:-1, :)
    reversed_weights = overlap_weights(3:1:-1)
    reversed_ownership = ownership(3:1:-1)
    call assemble_fci_boundary_patch_mortar( &
        background_trace_reversed, patch_trace_reversed, reversed_weights, &
        reversed_ownership, cross_mass_reversed, background_mass_reversed, &
        patch_mass_reversed, background_transfer_reversed, &
        patch_transfer_reversed, overlap_measure_reversed, status)
    call check_condition(status%code == 0, &
        "reversed interface orientation remains admissible")
    call check_condition(maxval(abs(cross_mass_reversed - cross_mass)) < tolerance, &
        "reversing quadrature orientation preserves the cross-mass")
    call check_condition(abs(overlap_measure_reversed - overlap_measure) < tolerance, &
        "reversing quadrature orientation preserves overlap measure")

    zero_weights = [1.0_dp, 0.0_dp, 1.0_dp]
    call assemble_fci_boundary_patch_mortar( &
        background_trace, patch_trace, zero_weights, ownership, cross_mass, &
        background_mass, patch_mass, background_from_patch, &
        patch_from_background, overlap_measure, status)
    call check_condition(status%code /= 0, &
        "FCI boundary-patch mortar rejects zero-measure quadrature")

    call assemble_fci_boundary_patch_mortar( &
        background_trace, patch_trace, overlap_weights, [1, 2, 1], cross_mass, &
        background_mass, patch_mass, background_from_patch, &
        patch_from_background, overlap_measure, status)
    call check_condition(status%code /= 0, &
        "FCI boundary-patch mortar rejects duplicate ownership")

    rank_deficient_trace(:, 1) = [0.5_dp, 0.5_dp, 0.5_dp]
    rank_deficient_trace(:, 2) = [0.5_dp, 0.5_dp, 0.5_dp]
    call assemble_fci_boundary_patch_mortar( &
        rank_deficient_trace, patch_trace, overlap_weights, ownership, cross_mass, &
        background_mass, patch_mass, background_from_patch, &
        patch_from_background, overlap_measure, status)
    call check_condition(status%code /= 0, &
        "FCI boundary-patch mortar rejects rank-deficient coupling")

    call check_summary("FCI boundary-patch mortar")
end program test_fci_boundary_patch_mortar
