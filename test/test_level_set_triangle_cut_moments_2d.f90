program test_level_set_triangle_cut_moments_2d
    use check, only: check_condition, check_summary
    use fortfem_core, only: evaluate_level_set_triangle_cut_moments_2d, &
        evaluate_level_set_triangle_cut_moments_2d_jvp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(2, 3) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
    real(dp), parameter :: level_values(3) = [1.0_dp, -1.0_dp, -1.0_dp]
    real(dp), parameter :: vertices_dot(2, 3) = reshape([ &
        0.03_dp, -0.02_dp, -0.01_dp, 0.04_dp, 0.02_dp, -0.01_dp], [2, 3])
    real(dp), parameter :: level_values_dot(3) = [0.07_dp, -0.03_dp, 0.05_dp]
    real(dp), parameter :: expected_centroid(2) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp]
    real(dp), parameter :: expected_second(2, 2) = reshape([ &
        1.0_dp/192.0_dp, 1.0_dp/384.0_dp, 1.0_dp/384.0_dp, &
        1.0_dp/192.0_dp], [2, 2])
    real(dp), parameter :: parent_second(2, 2) = reshape([ &
        1.0_dp/12.0_dp, 1.0_dp/24.0_dp, 1.0_dp/24.0_dp, &
        1.0_dp/12.0_dp], [2, 2])
    real(dp) :: positive_area, positive_centroid(2), positive_second(2, 2)
    real(dp) :: negative_area, negative_centroid(2), negative_second(2, 2)
    real(dp) :: interface_length, normal(2)
    real(dp) :: positive_area_dot, positive_centroid_dot(2), positive_second_dot(2, 2)
    real(dp) :: negative_area_dot, negative_centroid_dot(2), negative_second_dot(2, 2)
    real(dp) :: interface_length_dot, normal_dot(2)
    real(dp) :: plus_area, plus_centroid(2), plus_second(2, 2), plus_length, plus_normal(2)
    real(dp) :: minus_area, minus_centroid(2), minus_second(2, 2), minus_length, minus_normal(2)
    real(dp) :: fd_plus_area, fd_minus_area, fd_plus_centroid(2), fd_minus_centroid(2)
    real(dp) :: fd_plus_second(2, 2)
    real(dp) :: fd_second_dot(2, 2), fd_centroid_dot(2), fd_area_dot, fd_length_dot
    real(dp) :: fd_normal_dot(2), step
    integer :: status, plus_status, minus_status

    call evaluate_level_set_triangle_cut_moments_2d( &
        vertices, level_values, positive_area, positive_centroid, positive_second, &
        negative_area, negative_centroid, negative_second, interface_length, normal, &
        status)
    call check_condition(status == 0 .and. abs(positive_area - 1.0_dp/8.0_dp) < &
        1.0e-14_dp .and. maxval(abs(positive_centroid - expected_centroid)) < &
        1.0e-14_dp .and. maxval(abs(positive_second - expected_second)) < &
        1.0e-14_dp, "cut moments match the independent one-corner oracle")
    call check_condition(abs(positive_area + negative_area - 0.5_dp) < 1.0e-14_dp .and. &
        maxval(abs(positive_second + negative_second - parent_second)) < 1.0e-14_dp, &
        "cut moments conserve parent area and quadratic moments")

    call evaluate_level_set_triangle_cut_moments_2d( &
        vertices, [1.0_dp, 2.0_dp, 3.0_dp], positive_area, positive_centroid, &
        positive_second, negative_area, negative_centroid, negative_second, &
        interface_length, normal, status)
    call check_condition(status == 0 .and. abs(positive_area - 0.5_dp) < 1.0e-14_dp .and. &
        maxval(abs(positive_second - parent_second)) < 1.0e-14_dp .and. &
        abs(negative_area) < 1.0e-14_dp .and. maxval(abs(negative_second)) < &
        1.0e-14_dp, &
        "uncut moments reduce to the parent triangle")

    call evaluate_level_set_triangle_cut_moments_2d_jvp( &
        vertices, level_values, vertices_dot, level_values_dot, positive_area_dot, &
        positive_centroid_dot, positive_second_dot, negative_area_dot, &
        negative_centroid_dot, negative_second_dot, interface_length_dot, normal_dot, &
        status)
    step = 1.0e-6_dp
    call evaluate_level_set_triangle_cut_moments_2d( &
        vertices + step*vertices_dot, level_values + step*level_values_dot, &
        fd_plus_area, fd_plus_centroid, fd_plus_second, negative_area, &
        negative_centroid, negative_second, plus_length, plus_normal, plus_status)
    call evaluate_level_set_triangle_cut_moments_2d( &
        vertices - step*vertices_dot, level_values - step*level_values_dot, &
        minus_area, minus_centroid, minus_second, negative_area, negative_centroid, &
        negative_second, minus_length, minus_normal, minus_status)
    fd_area_dot = (fd_plus_area - minus_area)/(2.0_dp*step)
    fd_centroid_dot = (fd_plus_centroid - minus_centroid)/(2.0_dp*step)
    fd_second_dot = (fd_plus_second - minus_second)/(2.0_dp*step)
    fd_length_dot = (plus_length - minus_length)/(2.0_dp*step)
    fd_normal_dot = (plus_normal - minus_normal)/(2.0_dp*step)
    call check_condition(status == 0 .and. plus_status == 0 .and. minus_status == 0 .and. &
        abs(positive_area_dot - fd_area_dot) < 1.0e-7_dp .and. &
        maxval(abs(positive_centroid_dot - fd_centroid_dot)) < 1.0e-7_dp .and. &
        maxval(abs(positive_second_dot - fd_second_dot)) < 1.0e-7_dp .and. &
        abs(interface_length_dot - fd_length_dot) < 1.0e-7_dp .and. &
        maxval(abs(normal_dot - fd_normal_dot)) < 1.0e-7_dp, &
        "cut-moment JVP matches a fixed-topology central-difference oracle")

    call check_summary("2D level-set cut moments")
end program test_level_set_triangle_cut_moments_2d
