program test_fci_triangle_interpolation_map
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_fci_triangle_interpolation_map_2d, &
        build_fci_triangle_interpolation_map_2d_jvp, &
        build_fci_triangle_interpolation_map_2d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: vertices(2, 4) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, &
        0.0_dp, 1.0_dp], [2, 4])
    integer, parameter :: triangles(3, 2) = reshape([ &
        1, 2, 3, 1, 3, 4], [3, 2])
    real(dp), parameter :: target_points(2, 2) = reshape([ &
        0.40_dp, 0.20_dp, 0.75_dp, 0.90_dp], [2, 2])
    integer, parameter :: target_cells(2) = [1, 2]
    real(dp), parameter :: vertex_dot(2, 4) = reshape([ &
        0.01_dp, -0.02_dp, 0.03_dp, -0.04_dp, &
        0.05_dp, -0.06_dp, 0.07_dp, -0.08_dp], [2, 4])
    real(dp), parameter :: target_dot(2, 2) = reshape([ &
        0.04_dp, -0.03_dp, 0.02_dp, -0.01_dp], [2, 2])
    real(dp), parameter :: map_bar(2, 4) = reshape([ &
        0.7_dp, -0.3_dp, 0.2_dp, 0.1_dp, &
        0.6_dp, -0.4_dp, -0.1_dp, 0.5_dp], [2, 4])
    real(dp), parameter :: epsilon = 1.0e-6_dp
    real(dp) :: interpolation_map(2, 4), interpolation_map_dot(2, 4)
    real(dp) :: interpolation_map_plus(2, 4), interpolation_map_minus(2, 4)
    real(dp) :: vertices_bar(2, 4), target_points_bar(2, 2)
    real(dp) :: left_pairing, right_pairing
    real(dp), parameter :: affine_values(4) = [0.0_dp, 1.0_dp, 3.0_dp, 2.0_dp]
    real(dp) :: interpolated_values(2)
    integer, parameter :: bad_cells(2) = [1, 3]
    type(fortsparse_status_t) :: status

    call build_fci_triangle_interpolation_map_2d( &
        vertices, triangles, target_points, target_cells, interpolation_map, &
        status)
    call check_condition(status%code == 0, &
        "FCI triangle interpolation accepts fixed containing cells")
    interpolated_values = matmul(interpolation_map, affine_values)
    call check_condition(maxval(abs(interpolated_values - [0.80_dp, 2.55_dp])) < &
        1.0e-14_dp, "FCI triangle interpolation reproduces an affine field")
    call check_condition(maxval(abs(sum(interpolation_map, dim=2) - 1.0_dp)) < &
        1.0e-14_dp, "FCI triangle weights form a partition of unity")

    call build_fci_triangle_interpolation_map_2d_jvp( &
        vertices, triangles, target_points, target_cells, vertex_dot, target_dot, &
        interpolation_map_dot, status)
    call build_fci_triangle_interpolation_map_2d( &
        vertices + epsilon*vertex_dot, triangles, &
        target_points + epsilon*target_dot, target_cells, interpolation_map_plus, &
        status)
    call build_fci_triangle_interpolation_map_2d( &
        vertices - epsilon*vertex_dot, triangles, &
        target_points - epsilon*target_dot, target_cells, interpolation_map_minus, &
        status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(interpolation_map_dot - (interpolation_map_plus - &
        interpolation_map_minus)/(2.0_dp*epsilon))) < 1.0e-8_dp, &
        "FCI triangle JVP matches the central-difference oracle")

    call build_fci_triangle_interpolation_map_2d_vjp( &
        vertices, triangles, target_points, target_cells, map_bar, vertices_bar, &
        target_points_bar, status)
    left_pairing = sum(map_bar*interpolation_map_dot)
    right_pairing = sum(vertices_bar*vertex_dot) + &
        sum(target_points_bar*target_dot)
    call check_condition(status%code == 0 .and. &
        abs(left_pairing - right_pairing) < 1.0e-12_dp, &
        "FCI triangle VJP satisfies the real dot-product identity")

    call build_fci_triangle_interpolation_map_2d( &
        vertices, triangles, target_points, bad_cells, interpolation_map, status)
    call check_condition(status%code /= 0, &
        "FCI triangle interpolation rejects an invalid containing cell")
    call check_summary("FCI triangle interpolation map")
end program test_fci_triangle_interpolation_map
