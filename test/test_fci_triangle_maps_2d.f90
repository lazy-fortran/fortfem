program test_fci_triangle_maps_2d
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        build_fci_triangle_interpolation_maps_2d, &
        build_fci_triangle_interpolation_maps_2d_jvp, &
        build_fci_triangle_interpolation_maps_2d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: vertices(2, 4) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, &
        0.0_dp, 1.0_dp], [2, 4])
    integer, parameter :: triangles(3, 2) = reshape([ &
        1, 2, 3, 1, 3, 4], [3, 2])
    real(dp), parameter :: forward_points(2, 2, 2) = reshape([ &
        0.40_dp, 0.20_dp, 0.20_dp, 0.10_dp, &
        0.75_dp, 0.90_dp, 0.60_dp, 0.80_dp], [2, 2, 2])
    real(dp), parameter :: backward_points(2, 2, 2) = reshape([ &
        0.30_dp, 0.10_dp, 0.50_dp, 0.20_dp, &
        0.80_dp, 0.90_dp, 0.70_dp, 0.85_dp], [2, 2, 2])
    integer, parameter :: forward_cells(2, 2) = reshape([1, 1, 2, 2], [2, 2])
    integer, parameter :: backward_cells(2, 2) = reshape([1, 1, 2, 2], [2, 2])
    real(dp), parameter :: vertex_dot(2, 4) = reshape([ &
        0.01_dp, -0.02_dp, 0.03_dp, -0.04_dp, &
        0.05_dp, -0.06_dp, 0.07_dp, -0.08_dp], [2, 4])
    real(dp), parameter :: forward_points_dot(2, 2, 2) = reshape([ &
        0.04_dp, -0.03_dp, 0.02_dp, -0.01_dp, &
        0.03_dp, -0.02_dp, 0.01_dp, -0.04_dp], [2, 2, 2])
    real(dp), parameter :: backward_points_dot(2, 2, 2) = reshape([ &
        -0.03_dp, 0.02_dp, -0.01_dp, 0.04_dp, &
        -0.02_dp, 0.01_dp, -0.04_dp, 0.03_dp], [2, 2, 2])
    real(dp), parameter :: forward_map_bar(2, 4, 2) = reshape([ &
        0.7_dp, -0.3_dp, 0.2_dp, 0.1_dp, -0.1_dp, 0.6_dp, 0.4_dp, -0.2_dp, &
        0.6_dp, -0.4_dp, -0.1_dp, 0.5_dp, 0.2_dp, -0.5_dp, 0.4_dp, 0.8_dp], &
        [2, 4, 2])
    real(dp), parameter :: backward_map_bar(2, 4, 2) = reshape([ &
        -0.2_dp, 0.4_dp, -0.6_dp, 0.8_dp, -1.0_dp, 1.2_dp, 0.3_dp, -0.5_dp, &
        0.7_dp, -0.9_dp, 1.1_dp, -1.3_dp, 0.4_dp, -0.6_dp, 0.8_dp, -1.0_dp], &
        [2, 4, 2])
    real(dp), parameter :: epsilon = 1.0e-6_dp
    real(dp) :: forward_map(2, 4, 2), backward_map(2, 4, 2)
    real(dp) :: forward_map_dot(2, 4, 2), backward_map_dot(2, 4, 2)
    real(dp) :: forward_map_plus(2, 4, 2), backward_map_plus(2, 4, 2)
    real(dp) :: forward_map_minus(2, 4, 2), backward_map_minus(2, 4, 2)
    real(dp) :: vertices_bar(2, 4)
    real(dp) :: forward_points_bar(2, 2, 2), backward_points_bar(2, 2, 2)
    real(dp) :: left_pairing, right_pairing
    real(dp), parameter :: topology_points(2, 2, 2) = reshape([ &
        0.50_dp, 0.50_dp, 0.20_dp, 0.10_dp, &
        0.75_dp, 0.90_dp, 0.60_dp, 0.80_dp], [2, 2, 2])
    type(fortsparse_status_t) :: status

    call build_fci_triangle_interpolation_maps_2d( &
        vertices, triangles, forward_points, forward_cells, backward_points, &
        backward_cells, forward_map, backward_map, status)
    call check_condition(status%code == 0, &
        "FCI triangle map adapter accepts fixed per-segment cells")
    call check_condition(maxval(abs(sum(forward_map, dim=2) - 1.0_dp)) < &
        1.0e-14_dp .and. maxval(abs(sum(backward_map, dim=2) - 1.0_dp)) < &
        1.0e-14_dp, "FCI triangle map adapter preserves partition of unity")

    call build_fci_triangle_interpolation_maps_2d_jvp( &
        vertices, triangles, forward_points, forward_cells, backward_points, &
        backward_cells, vertex_dot, forward_points_dot, backward_points_dot, &
        forward_map_dot, backward_map_dot, status)
    call build_fci_triangle_interpolation_maps_2d( &
        vertices + epsilon*vertex_dot, triangles, &
        forward_points + epsilon*forward_points_dot, forward_cells, &
        backward_points + epsilon*backward_points_dot, backward_cells, &
        forward_map_plus, backward_map_plus, status)
    call build_fci_triangle_interpolation_maps_2d( &
        vertices - epsilon*vertex_dot, triangles, &
        forward_points - epsilon*forward_points_dot, forward_cells, &
        backward_points - epsilon*backward_points_dot, backward_cells, &
        forward_map_minus, backward_map_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(forward_map_dot - (forward_map_plus - forward_map_minus)/ &
        (2.0_dp*epsilon))) < 1.0e-8_dp .and. &
        maxval(abs(backward_map_dot - (backward_map_plus - backward_map_minus)/ &
        (2.0_dp*epsilon))) < 1.0e-8_dp, &
        "FCI triangle map adapter JVP matches the central-difference oracle")

    call build_fci_triangle_interpolation_maps_2d_vjp( &
        vertices, triangles, forward_points, forward_cells, backward_points, &
        backward_cells, forward_map_bar, backward_map_bar, vertices_bar, &
        forward_points_bar, backward_points_bar, status)
    left_pairing = sum(forward_map_bar*forward_map_dot) + &
        sum(backward_map_bar*backward_map_dot)
    right_pairing = sum(vertices_bar*vertex_dot) + &
        sum(forward_points_bar*forward_points_dot) + &
        sum(backward_points_bar*backward_points_dot)
    call check_condition(status%code == 0 .and. &
        abs(left_pairing - right_pairing) < 1.0e-12_dp, &
        "FCI triangle map adapter VJP satisfies the full dot-product identity")

    call build_fci_triangle_interpolation_maps_2d_jvp( &
        vertices, triangles, topology_points, forward_cells, backward_points, &
        backward_cells, vertex_dot, forward_points_dot, backward_points_dot, &
        forward_map_dot, backward_map_dot, status)
    call check_condition(status%code /= 0, &
        "FCI triangle map adapter JVP rejects a cell-boundary event")
    call check_summary("FCI triangle map adapter")
end program test_fci_triangle_maps_2d
