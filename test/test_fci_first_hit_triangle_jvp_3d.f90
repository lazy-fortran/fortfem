program test_fci_first_hit_triangle_jvp_3d
    use check, only: check_condition, check_summary
    use fortfem_fci_terminal_triangle_3d, only: &
        find_fci_first_hit_triangle_3d, find_fci_first_hit_triangle_3d_jvp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: surface_vertices(3, 6) = reshape([ &
        0.8_dp, 0.0_dp, 0.0_dp, 0.8_dp, 1.0_dp, 0.0_dp, &
        0.8_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], [3, 6])
    integer, parameter :: surface_triangles(3, 2) = reshape([ &
        1, 2, 3, 4, 5, 6], [3, 2])
    real(dp), parameter :: start_point(3) = [0.5_dp, 0.2_dp, 0.2_dp]
    real(dp), parameter :: finish_point(3) = [1.3_dp, 0.3_dp, 0.4_dp]
    real(dp), parameter :: start_point_dot(3) = [0.02_dp, -0.01_dp, 0.03_dp]
    real(dp), parameter :: finish_point_dot(3) = [-0.03_dp, 0.04_dp, -0.02_dp]
    real(dp), parameter :: surface_vertices_dot(3, 6) = reshape([ &
        0.01_dp, -0.02_dp, 0.03_dp, -0.01_dp, 0.02_dp, -0.04_dp, &
        0.02_dp, -0.01_dp, 0.03_dp, -0.02_dp, 0.01_dp, 0.04_dp, &
        -0.02_dp, 0.01_dp, -0.03_dp, 0.02_dp, -0.01_dp, 0.03_dp], [3, 6])
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: hit_point(3), hit_parameter, normal(3)
    real(dp) :: hit_point_dot(3), hit_parameter_dot, normal_dot(3)
    real(dp) :: hit_point_plus(3), hit_parameter_plus, normal_plus(3)
    real(dp) :: hit_point_minus(3), hit_parameter_minus, normal_minus(3)
    integer :: hit_triangle, status, status_plus, status_minus

    call find_fci_first_hit_triangle_3d( &
        start_point, finish_point, surface_vertices, surface_triangles, &
        hit_point, hit_parameter, hit_triangle, normal, status)
    call check_condition(status == 0 .and. hit_triangle == 1, &
        "FCI 3D first-hit JVP has a fixed terminal triangle")
    call find_fci_first_hit_triangle_3d_jvp( &
        start_point, finish_point, surface_vertices, surface_triangles, &
        start_point_dot, finish_point_dot, surface_vertices_dot, hit_point_dot, &
        hit_parameter_dot, normal_dot, status)
    call check_condition(status == 0, &
        "FCI 3D first-hit JVP accepts fixed terminal topology")

    call find_fci_first_hit_triangle_3d( &
        start_point + step*start_point_dot, finish_point + step*finish_point_dot, &
        surface_vertices + step*surface_vertices_dot, surface_triangles, &
        hit_point_plus, hit_parameter_plus, hit_triangle, normal_plus, status_plus)
    call find_fci_first_hit_triangle_3d( &
        start_point - step*start_point_dot, finish_point - step*finish_point_dot, &
        surface_vertices - step*surface_vertices_dot, surface_triangles, &
        hit_point_minus, hit_parameter_minus, hit_triangle, normal_minus, status_minus)
    call check_condition(status_plus == 0 .and. status_minus == 0 .and. &
        hit_triangle == 1, "FCI 3D first-hit central differences retain the triangle")
    call check_condition(maxval(abs(hit_point_dot - (hit_point_plus - &
        hit_point_minus)/(2.0_dp*step))) < 5.0e-9_dp .and. &
        abs(hit_parameter_dot - (hit_parameter_plus - hit_parameter_minus)/ &
        (2.0_dp*step)) < 5.0e-9_dp .and. maxval(abs(normal_dot - &
        (normal_plus - normal_minus)/(2.0_dp*step))) < 5.0e-9_dp, &
        "FCI 3D first-hit JVP matches independent central differences")

    call find_fci_first_hit_triangle_3d_jvp( &
        [0.5_dp, 1.2_dp, 1.2_dp], [1.3_dp, 1.3_dp, 1.4_dp], surface_vertices, &
        surface_triangles, start_point_dot, finish_point_dot, surface_vertices_dot, &
        hit_point_dot, hit_parameter_dot, normal_dot, status)
    call check_condition(status == 0 .and. maxval(abs(hit_point_dot)) == 0.0_dp .and. &
        hit_parameter_dot == 0.0_dp .and. maxval(abs(normal_dot)) == 0.0_dp, &
        "FCI 3D first-hit JVP handles a fixed no-hit trace")
    call check_summary("FCI first-hit terminal triangle JVP")
end program test_fci_first_hit_triangle_jvp_3d
