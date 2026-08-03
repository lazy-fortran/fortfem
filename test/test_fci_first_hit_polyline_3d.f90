program test_fci_first_hit_polyline_3d
    use check, only: check_condition, check_summary
    use fortfem_fci_terminal_polyline_3d, only: find_fci_first_hit_polyline_3d
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
    real(dp) :: hit_point(3), normal(3), connection_length
    real(dp) :: expected_length, first_length, second_length
    integer :: hit_segment, hit_triangle, status

    first_length = norm2(polyline(:, 2) - polyline(:, 1))
    second_length = norm2(polyline(:, 3) - polyline(:, 2))
    expected_length = first_length + second_length/6.0_dp
    call find_fci_first_hit_polyline_3d( &
        polyline, surface_vertices, surface_triangles, hit_point, connection_length, &
        hit_segment, hit_triangle, normal, status)
    call check_condition(status == 0 .and. hit_segment == 2 .and. hit_triangle == 1, &
        "FCI polyline search returns the first hit segment and triangle")
    call check_condition(maxval(abs(hit_point - &
        [0.8_dp, 0.233333333333333_dp, 0.266666666666667_dp])) < 1.0e-13_dp .and. &
        abs(connection_length - expected_length) < 1.0e-13_dp .and. &
        maxval(abs(normal - [1.0_dp, 0.0_dp, 0.0_dp])) < 1.0e-14_dp, &
        "FCI polyline search returns hit point, connection length, and normal")

    call find_fci_first_hit_polyline_3d( &
        no_hit_polyline, surface_vertices, surface_triangles, hit_point, &
        connection_length, hit_segment, hit_triangle, normal, status)
    call check_condition(status == 0 .and. hit_segment == 0 .and. &
        hit_triangle == 0 .and. &
        maxval(abs(hit_point - no_hit_polyline(:, 3))) < 1.0e-14_dp .and. &
        abs(connection_length - &
        (norm2(no_hit_polyline(:, 2) - no_hit_polyline(:, 1)) + &
        norm2(no_hit_polyline(:, 3) - no_hit_polyline(:, 2)))) < 1.0e-14_dp .and. &
        maxval(abs(normal)) == 0.0_dp, &
        "FCI polyline search returns endpoint and length for a valid no-hit trace")

    call find_fci_first_hit_polyline_3d( &
        polyline, surface_vertices, reshape([1, 1, 3, 4, 5, 6], [3, 2]), &
        hit_point, connection_length, hit_segment, hit_triangle, normal, status)
    call check_condition(status /= 0, "FCI polyline search rejects a degenerate facet")
    call check_summary("FCI first-hit terminal polyline")
end program test_fci_first_hit_polyline_3d
