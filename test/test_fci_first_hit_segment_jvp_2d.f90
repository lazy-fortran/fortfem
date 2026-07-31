program test_fci_first_hit_segment_jvp_2d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        find_fci_first_hit_segment_2d, find_fci_first_hit_segment_2d_jvp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: wall_vertices(2, 6) = reshape([ &
        0.8_dp, 0.0_dp, 0.8_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
        1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 6])
    integer, parameter :: wall_segments(2, 3) = reshape([ &
        1, 2, 3, 4, 5, 6], [2, 3])
    real(dp), parameter :: start_point(2) = [0.5_dp, 0.5_dp]
    real(dp), parameter :: finish_point(2) = [1.2_dp, 0.75_dp]
    real(dp), parameter :: start_point_dot(2) = [0.02_dp, -0.01_dp]
    real(dp), parameter :: finish_point_dot(2) = [-0.03_dp, 0.04_dp]
    real(dp), parameter :: wall_vertices_dot(2, 6) = reshape([ &
        0.01_dp, -0.02_dp, 0.0_dp, 0.03_dp, -0.01_dp, 0.02_dp, &
        0.02_dp, -0.01_dp, 0.03_dp, -0.02_dp, 0.01_dp, -0.04_dp], [2, 6])
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: hit_point(2), hit_parameter, normal(2)
    real(dp) :: hit_point_dot(2), hit_parameter_dot, normal_dot(2)
    real(dp) :: hit_point_plus(2), hit_parameter_plus, normal_plus(2)
    real(dp) :: hit_point_minus(2), hit_parameter_minus, normal_minus(2)
    integer :: hit_segment, status, status_plus, status_minus

    call find_fci_first_hit_segment_2d( &
        start_point, finish_point, wall_vertices, wall_segments, hit_point, &
        hit_parameter, hit_segment, normal, status)
    call check_condition(status == 0 .and. hit_segment == 1, &
        "FCI first-hit JVP has a fixed terminal facet")
    call find_fci_first_hit_segment_2d_jvp( &
        start_point, finish_point, wall_vertices, wall_segments, start_point_dot, &
        finish_point_dot, wall_vertices_dot, hit_point_dot, hit_parameter_dot, &
        normal_dot, status)
    call check_condition(status == 0, &
        "FCI first-hit JVP accepts fixed terminal topology")

    call find_fci_first_hit_segment_2d( &
        start_point + step*start_point_dot, finish_point + step*finish_point_dot, &
        wall_vertices + step*wall_vertices_dot, wall_segments, hit_point_plus, &
        hit_parameter_plus, hit_segment, normal_plus, status_plus)
    call find_fci_first_hit_segment_2d( &
        start_point - step*start_point_dot, finish_point - step*finish_point_dot, &
        wall_vertices - step*wall_vertices_dot, wall_segments, hit_point_minus, &
        hit_parameter_minus, hit_segment, normal_minus, status_minus)
    call check_condition(status_plus == 0 .and. status_minus == 0 .and. &
        hit_segment == 1, "FCI first-hit central differences retain the facet")
    call check_condition(maxval(abs(hit_point_dot - (hit_point_plus - &
        hit_point_minus)/(2.0_dp*step))) < 5.0e-9_dp .and. &
        abs(hit_parameter_dot - (hit_parameter_plus - hit_parameter_minus)/ &
        (2.0_dp*step)) < 5.0e-9_dp .and. maxval(abs(normal_dot - &
        (normal_plus - normal_minus)/(2.0_dp*step))) < 5.0e-9_dp, &
        "FCI first-hit JVP matches independent central differences")

    call find_fci_first_hit_segment_jvp_no_hit( &
        hit_point_dot, hit_parameter_dot, normal_dot, status)
    call check_condition(status == 0 .and. maxval(abs(hit_point_dot)) == 0.0_dp .and. &
        hit_parameter_dot == 0.0_dp .and. maxval(abs(normal_dot)) == 0.0_dp, &
        "FCI first-hit JVP handles a fixed no-hit trace")
    call check_summary("FCI first-hit terminal segment JVP")

contains

    subroutine find_fci_first_hit_segment_jvp_no_hit( &
            hit_point_dot, hit_parameter_dot, normal_dot, status)
        real(dp), intent(out) :: hit_point_dot(2), hit_parameter_dot, normal_dot(2)
        integer, intent(out) :: status

        call find_fci_first_hit_segment_2d_jvp( &
            [0.5_dp, 0.5_dp], [0.7_dp, 0.7_dp], wall_vertices, wall_segments, &
            [0.0_dp, 0.01_dp], [0.02_dp, -0.01_dp], wall_vertices_dot, &
            hit_point_dot, hit_parameter_dot, normal_dot, status)
    end subroutine find_fci_first_hit_segment_jvp_no_hit

end program test_fci_first_hit_segment_jvp_2d
