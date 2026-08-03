program test_fci_first_hit_triangle_3d
    use check, only: check_condition, check_summary
    use fortfem_fci_terminal_triangle_3d, only: find_fci_first_hit_triangle_3d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: surface_vertices(3, 6) = reshape([ &
        0.8_dp, 0.0_dp, 0.0_dp, 0.8_dp, 1.0_dp, 0.0_dp, &
        0.8_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], [3, 6])
    integer, parameter :: surface_triangles(3, 2) = reshape([ &
        1, 2, 3, 4, 5, 6], [3, 2])
    real(dp) :: hit_point(3), normal(3), hit_parameter
    integer :: hit_triangle, status

    call find_fci_first_hit_triangle_3d( &
        [0.5_dp, 0.2_dp, 0.2_dp], [1.3_dp, 0.3_dp, 0.4_dp], &
        surface_vertices, surface_triangles, hit_point, hit_parameter, &
        hit_triangle, normal, status)
    call check_condition(status == 0 .and. hit_triangle == 1, &
        "FCI 3D first-hit search returns the nearest triangle")
    call check_condition(abs(hit_parameter - 0.375_dp) < 1.0e-14_dp .and. &
        maxval(abs(hit_point - [0.8_dp, 0.2375_dp, 0.275_dp])) < 1.0e-14_dp .and. &
        maxval(abs(normal - [1.0_dp, 0.0_dp, 0.0_dp])) < 1.0e-14_dp, &
        "FCI 3D first-hit search returns exact point and oriented normal")

    call find_fci_first_hit_triangle_3d( &
        [0.5_dp, 1.2_dp, 1.2_dp], [1.3_dp, 1.3_dp, 1.4_dp], &
        surface_vertices, surface_triangles, hit_point, hit_parameter, &
        hit_triangle, normal, status)
    call check_condition(status == 0 .and. hit_triangle == 0 .and. &
        hit_parameter == 0.0_dp .and. maxval(abs(hit_point)) == 0.0_dp .and. &
        maxval(abs(normal)) == 0.0_dp, &
        "FCI 3D first-hit search reports a valid trace with no terminal hit")

    call find_fci_first_hit_triangle_3d( &
        [0.5_dp, 0.2_dp, 0.2_dp], [1.3_dp, 0.3_dp, 0.4_dp], &
        surface_vertices, reshape([2, 1, 3, 4, 5, 6], [3, 2]), hit_point, &
        hit_parameter, hit_triangle, normal, status)
    call check_condition(status == 0 .and. hit_triangle == 1 .and. &
        maxval(abs(normal - [-1.0_dp, 0.0_dp, 0.0_dp])) < 1.0e-14_dp, &
        "FCI 3D first-hit search preserves triangle orientation")

    call find_fci_first_hit_triangle_3d( &
        [0.5_dp, 0.2_dp, 0.2_dp], [1.3_dp, 0.3_dp, 0.4_dp], &
        surface_vertices, reshape([1, 1, 3, 4, 5, 6], [3, 2]), hit_point, &
        hit_parameter, hit_triangle, normal, status)
    call check_condition(status /= 0, &
        "FCI 3D first-hit search rejects a degenerate triangle")
    call check_summary("FCI first-hit terminal triangle")
end program test_fci_first_hit_triangle_3d
