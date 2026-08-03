program test_fci_first_hit_segment_2d
    use check, only: check_condition, check_summary
    use fortfem_fci_terminal_segment_2d, only: find_fci_first_hit_segment_2d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: wall_vertices(2, 6) = reshape([ &
        0.8_dp, 0.0_dp, 0.8_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
        1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 6])
    integer, parameter :: wall_segments(2, 3) = reshape([ &
        1, 2, 3, 4, 5, 6], [2, 3])
    real(dp) :: hit_point(2), normal(2), hit_parameter
    integer :: hit_segment, status

    call find_fci_first_hit_segment_2d( &
        [0.5_dp, 0.5_dp], [1.2_dp, 0.75_dp], wall_vertices, wall_segments, &
        hit_point, hit_parameter, hit_segment, normal, status)
    call check_condition(status == 0 .and. hit_segment == 1, &
        "FCI first-hit search returns the nearest oriented wall segment")
    call check_condition(abs(hit_parameter - 3.0_dp/7.0_dp) < 1.0e-14_dp .and. &
        maxval(abs(hit_point - [0.8_dp, 17.0_dp/28.0_dp])) < 1.0e-14_dp .and. &
        maxval(abs(normal - [1.0_dp, 0.0_dp])) < 1.0e-14_dp, &
        "FCI first-hit search returns the exact intersection and right normal")

    call find_fci_first_hit_segment_2d( &
        [0.5_dp, 0.5_dp], [0.7_dp, 0.7_dp], wall_vertices, wall_segments, &
        hit_point, hit_parameter, hit_segment, normal, status)
    call check_condition(status == 0 .and. hit_segment == 0 .and. &
        hit_parameter == 0.0_dp .and. maxval(abs(hit_point)) == 0.0_dp .and. &
        maxval(abs(normal)) == 0.0_dp, &
        "FCI first-hit search reports a valid trace with no terminal hit")

    call find_fci_first_hit_segment_2d( &
        [0.5_dp, 0.5_dp], [1.2_dp, 0.75_dp], wall_vertices, &
        reshape([2, 1, 3, 4, 5, 6], [2, 3]), hit_point, hit_parameter, &
        hit_segment, normal, status)
    call check_condition(status == 0 .and. hit_segment == 1 .and. &
        maxval(abs(normal - [-1.0_dp, 0.0_dp])) < 1.0e-14_dp, &
        "FCI first-hit search preserves facet-local orientation")

    call find_fci_first_hit_segment_2d( &
        [0.5_dp, 0.5_dp], [1.2_dp, 0.75_dp], wall_vertices, &
        reshape([1, 1, 3, 4, 5, 6], [2, 3]), hit_point, &
        hit_parameter, hit_segment, normal, status)
    call check_condition(status /= 0, &
        "FCI first-hit search rejects a degenerate wall facet")
    call check_summary("FCI first-hit terminal segment")
end program test_fci_first_hit_segment_2d
