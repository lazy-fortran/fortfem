program test_fci_first_hit_polyline_jvp_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        find_fci_first_hit_polyline_3d, find_fci_first_hit_polyline_3d_jvp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: surface_vertices(3, 6) = reshape([ &
        0.8_dp, 0.0_dp, 0.0_dp, 0.8_dp, 1.0_dp, 0.0_dp, &
        0.8_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], [3, 6])
    integer, parameter :: surface_triangles(3, 2) = reshape([ &
        1, 2, 3, 4, 5, 6], [3, 2])
    real(dp), parameter :: polyline(3, 3) = reshape([ &
        0.5_dp, 0.2_dp, 0.2_dp, 0.7_dp, 0.22_dp, 0.24_dp, &
        1.3_dp, 0.3_dp, 0.4_dp], [3, 3])
    real(dp), parameter :: no_hit_polyline(3, 3) = reshape([ &
        0.5_dp, 1.2_dp, 1.2_dp, 0.7_dp, 1.22_dp, 1.24_dp, &
        1.3_dp, 1.3_dp, 1.4_dp], [3, 3])
    real(dp), parameter :: polyline_dot(3, 3) = reshape([ &
        0.02_dp, -0.01_dp, 0.03_dp, -0.03_dp, 0.04_dp, -0.02_dp, &
        0.01_dp, -0.02_dp, 0.04_dp], [3, 3])
    real(dp), parameter :: surface_vertices_dot(3, 6) = reshape([ &
        0.01_dp, -0.02_dp, 0.03_dp, -0.01_dp, 0.02_dp, -0.04_dp, &
        0.02_dp, -0.01_dp, 0.03_dp, -0.02_dp, 0.01_dp, 0.04_dp, &
        -0.02_dp, 0.01_dp, -0.03_dp, 0.02_dp, -0.01_dp, 0.03_dp], [3, 6])
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: hit_point(3), normal(3), connection_length
    real(dp) :: hit_point_dot(3), normal_dot(3), connection_length_dot
    real(dp) :: hit_point_plus(3), normal_plus(3), connection_length_plus
    real(dp) :: hit_point_minus(3), normal_minus(3), connection_length_minus
    integer :: hit_segment, hit_triangle, status, status_plus, status_minus

    call find_fci_first_hit_polyline_3d( &
        polyline, surface_vertices, surface_triangles, hit_point, connection_length, &
        hit_segment, hit_triangle, normal, status)
    call check_condition(status == 0 .and. hit_segment == 2 .and. hit_triangle == 1, &
        "FCI polyline JVP has a fixed terminal event")
    call find_fci_first_hit_polyline_3d_jvp( &
        polyline, surface_vertices, surface_triangles, polyline_dot, &
        surface_vertices_dot, hit_point_dot, connection_length_dot, normal_dot, status)
    call check_condition(status == 0, "FCI polyline JVP accepts fixed topology")

    call find_fci_first_hit_polyline_3d( &
        polyline + step*polyline_dot, surface_vertices + step*surface_vertices_dot, &
        surface_triangles, hit_point_plus, connection_length_plus, hit_segment, &
        hit_triangle, normal_plus, status_plus)
    call find_fci_first_hit_polyline_3d( &
        polyline - step*polyline_dot, surface_vertices - step*surface_vertices_dot, &
        surface_triangles, hit_point_minus, connection_length_minus, hit_segment, &
        hit_triangle, normal_minus, status_minus)
    call check_condition(status_plus == 0 .and. status_minus == 0 .and. &
        hit_segment == 2 .and. hit_triangle == 1, &
        "FCI polyline central differences retain the terminal event")
    call check_condition(maxval(abs(hit_point_dot - (hit_point_plus - &
        hit_point_minus)/ &
        (2.0_dp*step))) < 5.0e-9_dp .and. &
        abs(connection_length_dot - (connection_length_plus - &
        connection_length_minus)/ &
        (2.0_dp*step)) < 5.0e-9_dp .and. maxval(abs(normal_dot - &
        (normal_plus - normal_minus)/(2.0_dp*step))) < 5.0e-9_dp, &
        "FCI polyline JVP matches independent central differences")

    call find_fci_first_hit_polyline_3d_jvp( &
        no_hit_polyline, surface_vertices, surface_triangles, polyline_dot, &
        surface_vertices_dot, hit_point_dot, &
        connection_length_dot, normal_dot, status)
    call check_condition(status == 0, "FCI polyline JVP handles a no-hit path")
    call check_summary("FCI first-hit terminal polyline JVP")

end program test_fci_first_hit_polyline_jvp_3d
