program test_fci_parallel_operator
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_fci_parallel_gradient_csc, &
        assemble_fci_parallel_support_divergence_csc, &
        apply_fci_parallel_gradient, apply_fci_parallel_gradient_jvp, &
        apply_fci_parallel_gradient_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    real(dp), parameter :: forward_map(2, 2, 2) = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2, 2])
    real(dp), parameter :: backward_map(2, 2, 2) = forward_map
    real(dp), parameter :: line_lengths(2, 2) = reshape([ &
        1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [2, 2])
    real(dp), parameter :: canonical_volumes(6) = [ &
        1.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp]
    real(dp), parameter :: staggered_volumes(4) = [ &
        1.5_dp, 0.5_dp, 2.0_dp, 0.75_dp]
    real(dp), parameter :: field(6) = [ &
        0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: constant_field(6) = 7.0_dp
    real(dp), parameter :: flux(4) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    real(dp), parameter :: expected_divergence(6) = [ &
        1.5_dp, 0.5_dp, 4.5_dp, 2.0_dp, -3.0_dp, -3.0_dp]
    real(dp), allocatable :: gradient_field(:), constant_gradient(:)
    real(dp), allocatable :: divergence_flux(:)
    real(dp), allocatable :: gradient_dot(:), gradient_plus(:), gradient_minus(:)
    real(dp), allocatable :: forward_map_dot(:, :, :), backward_map_dot(:, :, :)
    real(dp), allocatable :: line_lengths_dot(:, :), field_dot(:)
    real(dp), allocatable :: bad_field_dot(:)
    real(dp), allocatable :: forward_map_bar(:, :, :), backward_map_bar(:, :, :)
    real(dp), allocatable :: line_lengths_bar(:, :), field_bar(:)
    real(dp), parameter :: finite_difference_step = 1.0e-7_dp
    real(dp), parameter :: output_bar(4) = [0.7_dp, -0.2_dp, 0.3_dp, 0.9_dp]
    real(dp), parameter :: map_direction(2, 2, 2) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp, &
        -0.5_dp, 0.6_dp, -0.7_dp, 0.8_dp], [2, 2, 2])
    real(dp), parameter :: reverse_map_direction(2, 2, 2) = reshape([ &
        -0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp, &
        0.7_dp, -0.8_dp, 0.9_dp, -1.0_dp], [2, 2, 2])
    real(dp), parameter :: length_direction(2, 2) = reshape([ &
        0.2_dp, -0.1_dp, 0.3_dp, -0.4_dp], [2, 2])
    real(dp), parameter :: field_direction(6) = [ &
        -0.2_dp, 0.4_dp, -0.6_dp, 0.8_dp, -1.0_dp, 1.2_dp]
    real(dp) :: vjp_left, vjp_right
    real(dp) :: weighted_pairing, adjoint_pairing
    type(csc_t) :: gradient, divergence
    type(fortsparse_status_t) :: status

    call assemble_fci_parallel_gradient_csc( &
        forward_map, backward_map, line_lengths, gradient, status)
    call check_condition(status%code == 0, &
        "FCI gradient assembly accepts two mapped plane interfaces")
    call check_condition(gradient%nrow == 4 .and. gradient%ncol == 6, &
        "FCI gradient has staggered rows and canonical columns")

    allocate(gradient_field(gradient%nrow), constant_gradient(gradient%nrow))
    gradient_field = csc_matvec(gradient, field)
    constant_gradient = csc_matvec(gradient, constant_field)
    call check_condition(maxval(abs(gradient_field - 1.0_dp)) < 1.0e-13_dp, &
        "FCI gradient reproduces a linear field along a straight line")
    call check_condition(maxval(abs(constant_gradient)) < 1.0e-13_dp, &
        "FCI gradient annihilates a constant field")

    allocate( &
        forward_map_dot(2, 2, 2), backward_map_dot(2, 2, 2), &
        line_lengths_dot(2, 2), field_dot(6), gradient_dot(4), &
        gradient_plus(4), gradient_minus(4))
    forward_map_dot = map_direction
    backward_map_dot = reverse_map_direction
    line_lengths_dot = length_direction
    field_dot = field_direction
    call apply_fci_parallel_gradient( &
        forward_map, backward_map, line_lengths, field, gradient_field, status)
    call apply_fci_parallel_gradient_jvp( &
        forward_map, backward_map, line_lengths, field, forward_map_dot, &
        backward_map_dot, line_lengths_dot, field_dot, gradient_dot, status)
    call apply_fci_parallel_gradient( &
        forward_map + finite_difference_step*forward_map_dot, &
        backward_map + finite_difference_step*backward_map_dot, &
        line_lengths + finite_difference_step*line_lengths_dot, &
        field + finite_difference_step*field_dot, gradient_plus, status)
    call apply_fci_parallel_gradient( &
        forward_map - finite_difference_step*forward_map_dot, &
        backward_map - finite_difference_step*backward_map_dot, &
        line_lengths - finite_difference_step*line_lengths_dot, &
        field - finite_difference_step*field_dot, gradient_minus, status)
    call check_condition(maxval(abs( &
        (gradient_plus - gradient_minus)/(2.0_dp*finite_difference_step) - &
        gradient_dot)) < 2.0e-8_dp, &
        "FCI gradient JVP matches an independent central difference")

    allocate(bad_field_dot(5))
    call apply_fci_parallel_gradient_jvp( &
        forward_map, backward_map, line_lengths, field, forward_map_dot, &
        backward_map_dot, line_lengths_dot, bad_field_dot, gradient_dot, status)
    call check_condition(status%code /= 0, &
        "FCI gradient JVP rejects an incompatible field tangent")

    allocate( &
        forward_map_bar(2, 2, 2), backward_map_bar(2, 2, 2), &
        line_lengths_bar(2, 2), field_bar(6))
    call apply_fci_parallel_gradient_vjp( &
        forward_map, backward_map, line_lengths, field, output_bar, &
        forward_map_bar, backward_map_bar, line_lengths_bar, field_bar, status)
    vjp_left = dot_product(output_bar, gradient_dot)
    vjp_right = sum(forward_map_bar*forward_map_dot) + &
        sum(backward_map_bar*backward_map_dot) + &
        sum(line_lengths_bar*line_lengths_dot) + dot_product(field_bar, field_dot)
    call check_condition(abs(vjp_left - vjp_right) < 2.0e-13_dp, &
        "FCI gradient VJP satisfies the real adjoint identity")

    call assemble_fci_parallel_support_divergence_csc( &
        gradient, canonical_volumes, staggered_volumes, divergence, status)
    call check_condition(status%code == 0 .and. divergence%nrow == 6 .and. &
        divergence%ncol == 4, &
        "FCI support divergence has the weighted adjoint shape")

    allocate(divergence_flux(divergence%nrow))
    divergence_flux = csc_matvec(divergence, flux)
    call check_condition(maxval(abs(divergence_flux - expected_divergence)) < &
        1.0e-13_dp, &
        "FCI support divergence matches the independent flux-balance oracle")

    weighted_pairing = dot_product( &
        canonical_volumes*divergence_flux, field)
    adjoint_pairing = -dot_product( &
        staggered_volumes*gradient_field, flux)
    call check_condition(abs(weighted_pairing - adjoint_pairing) < 1.0e-13_dp, &
        "FCI support divergence is the negative volume-weighted adjoint")

    call check_summary("PARALLAX-aligned FCI support operator")
end program test_fci_parallel_operator
