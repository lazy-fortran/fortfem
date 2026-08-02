program test_fci_reproducible_map
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_fci_parallel_gradient_csc, &
        assemble_fci_parallel_support_divergence_csc, &
        build_fci_bilinear_interpolation_maps_2d
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: nx = 3, ny = 2, plane_count = 4
    integer, parameter :: staggered_count = 2, segment_count = plane_count - 1
    integer, parameter :: plane_dof_count = nx*ny
    real(dp), parameter :: source_x(nx) = [0.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: source_y(ny) = [0.0_dp, 1.0_dp]
    real(dp), parameter :: forward_x(staggered_count, segment_count) = reshape([ &
        0.4_dp, 1.3_dp, 0.6_dp, 1.5_dp, 0.8_dp, 1.7_dp], &
        [staggered_count, segment_count])
    real(dp), parameter :: forward_y(staggered_count, segment_count) = reshape([ &
        0.2_dp, 0.7_dp, 0.3_dp, 0.6_dp, 0.4_dp, 0.5_dp], &
        [staggered_count, segment_count])
    real(dp), parameter :: backward_x(staggered_count, segment_count) = reshape([ &
        0.1_dp, 1.0_dp, 0.3_dp, 1.2_dp, 0.5_dp, 1.4_dp], &
        [staggered_count, segment_count])
    real(dp), parameter :: backward_y(staggered_count, segment_count) = reshape([ &
        0.1_dp, 0.8_dp, 0.2_dp, 0.7_dp, 0.3_dp, 0.6_dp], &
        [staggered_count, segment_count])
    real(dp), parameter :: line_lengths(staggered_count, segment_count) = reshape([ &
        0.5_dp, 0.75_dp, 0.6_dp, 0.8_dp, 0.7_dp, 0.9_dp], &
        [staggered_count, segment_count])
    real(dp), parameter :: canonical_volumes(plane_dof_count*plane_count) = [ &
        1.0_dp, 1.1_dp, 1.2_dp, 1.3_dp, 1.4_dp, 1.5_dp, &
        1.6_dp, 1.7_dp, 1.8_dp, 1.9_dp, 2.0_dp, 2.1_dp, &
        2.2_dp, 2.3_dp, 2.4_dp, 2.5_dp, 2.6_dp, 2.7_dp, &
        2.8_dp, 2.9_dp, 3.0_dp, 3.1_dp, 3.2_dp, 3.3_dp]
    real(dp), parameter :: staggered_volumes(staggered_count*segment_count) = &
        [1.2_dp, 1.4_dp, 1.6_dp, 1.8_dp, 2.0_dp, 2.2_dp]

    real(dp) :: forward_map(staggered_count, plane_dof_count, segment_count)
    real(dp) :: backward_map(staggered_count, plane_dof_count, segment_count)
    real(dp) :: forward_map_repeat(staggered_count, plane_dof_count, segment_count)
    real(dp) :: backward_map_repeat(staggered_count, plane_dof_count, segment_count)
    real(dp) :: field(plane_dof_count*plane_count)
    real(dp) :: gradient_field(staggered_count*segment_count)
    real(dp) :: expected_gradient(staggered_count*segment_count)
    real(dp) :: flux(staggered_count*segment_count) = [ &
        0.4_dp, -0.7_dp, 0.8_dp, -0.2_dp, 1.1_dp, -0.5_dp]
    real(dp) :: divergence_flux(plane_dof_count*plane_count)
    real(dp) :: weighted_pairing, adjoint_pairing
    type(csc_t) :: gradient, divergence
    type(fortsparse_status_t) :: status
    integer :: plane, i, j, node, segment, sample, row

    call build_fci_bilinear_interpolation_maps_2d( &
        source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
        forward_map, backward_map, status)
    call check_condition(status%code == 0, &
        "reproducible FCI map accepts deterministic traced endpoints")
    call check_condition(maxval(abs(sum(forward_map, dim=2) - 1.0_dp)) < 1.0e-14_dp .and. &
        maxval(abs(sum(backward_map, dim=2) - 1.0_dp)) < 1.0e-14_dp, &
        "reproducible FCI maps preserve partition of unity")
    call build_fci_bilinear_interpolation_maps_2d( &
        source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
        forward_map_repeat, backward_map_repeat, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(forward_map - forward_map_repeat)) == 0.0_dp .and. &
        maxval(abs(backward_map - backward_map_repeat)) == 0.0_dp, &
        "rebuilding the fixed FCI map is bitwise reproducible")

    do plane = 1, plane_count
        do j = 1, ny
            do i = 1, nx
                node = (plane - 1)*plane_dof_count + (j - 1)*nx + i
                field(node) = source_x(i) + 2.0_dp*source_y(j) + &
                    4.0_dp*real(plane - 1, dp)
            end do
        end do
    end do
    do segment = 1, segment_count
        do sample = 1, staggered_count
            row = sample + (segment - 1)*staggered_count
            expected_gradient(row) = ( &
                forward_x(sample, segment) + 2.0_dp*forward_y(sample, segment) + &
                4.0_dp*real(segment, dp) - backward_x(sample, segment) - &
                2.0_dp*backward_y(sample, segment) - &
                4.0_dp*real(segment - 1, dp))/line_lengths(sample, segment)
        end do
    end do

    call assemble_fci_parallel_gradient_csc( &
        forward_map, backward_map, line_lengths, gradient, status)
    call check_condition(status%code == 0, &
        "reproducible FCI map assembles the sparse parallel gradient")
    gradient_field = csc_matvec(gradient, field)
    call check_condition(maxval(abs(gradient_field - expected_gradient)) < 1.0e-13_dp, &
        "FCI gradient matches the independent affine field-line oracle")

    call assemble_fci_parallel_support_divergence_csc( &
        gradient, canonical_volumes, staggered_volumes, divergence, status)
    call check_condition(status%code == 0, &
        "reproducible FCI map assembles the weighted support divergence")
    divergence_flux = csc_matvec(divergence, flux)
    weighted_pairing = dot_product(canonical_volumes*divergence_flux, field)
    adjoint_pairing = -dot_product(staggered_volumes*gradient_field, flux)
    call check_condition(abs(weighted_pairing - adjoint_pairing) < 1.0e-13_dp, &
        "reproducible FCI map preserves the weighted support adjoint identity")

    call check_summary("reproducible FCI map")
end program test_fci_reproducible_map
