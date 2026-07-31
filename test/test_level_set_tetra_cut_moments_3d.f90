program test_level_set_tetra_cut_moments_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_level_set_tetra_cut_moments_3d, &
        evaluate_level_set_tetra_cut_moments_3d_jvp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(3, 4) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
    real(dp), parameter :: level_values(4) = [1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp]
    real(dp), parameter :: vertices_dot(3, 4) = reshape([ &
        0.02_dp, -0.01_dp, 0.03_dp, -0.02_dp, 0.04_dp, -0.01_dp, &
        0.01_dp, -0.03_dp, 0.02_dp, 0.03_dp, -0.02_dp, 0.01_dp], [3, 4])
    real(dp), parameter :: level_values_dot(4) = [0.06_dp, -0.02_dp, 0.04_dp, -0.03_dp]
    real(dp), parameter :: expected_centroid(3) = [1.0_dp/8.0_dp, 1.0_dp/8.0_dp, &
        1.0_dp/8.0_dp]
    real(dp), parameter :: expected_second(3, 3) = reshape([ &
        1.0_dp/1920.0_dp, 1.0_dp/3840.0_dp, 1.0_dp/3840.0_dp, &
        1.0_dp/3840.0_dp, 1.0_dp/1920.0_dp, 1.0_dp/3840.0_dp, &
        1.0_dp/3840.0_dp, 1.0_dp/3840.0_dp, 1.0_dp/1920.0_dp], [3, 3])
    real(dp), parameter :: parent_second(3, 3) = reshape([ &
        1.0_dp/60.0_dp, 1.0_dp/120.0_dp, 1.0_dp/120.0_dp, &
        1.0_dp/120.0_dp, 1.0_dp/60.0_dp, 1.0_dp/120.0_dp, &
        1.0_dp/120.0_dp, 1.0_dp/120.0_dp, 1.0_dp/60.0_dp], [3, 3])
    real(dp) :: positive_volume, positive_centroid(3), positive_second(3, 3)
    real(dp) :: negative_volume, negative_centroid(3), negative_second(3, 3)
    real(dp) :: interface_area, normal(3)
    real(dp) :: positive_volume_dot, positive_centroid_dot(3)
    real(dp) :: positive_second_dot(3, 3), negative_volume_dot
    real(dp) :: negative_centroid_dot(3), negative_second_dot(3, 3)
    real(dp) :: interface_area_dot, normal_dot(3)
    real(dp) :: plus_volume, plus_centroid(3), plus_second(3, 3)
    real(dp) :: minus_volume, minus_centroid(3), minus_second(3, 3)
    real(dp) :: plus_negative_volume, plus_negative_centroid(3)
    real(dp) :: plus_negative_second(3, 3), plus_area, plus_normal(3)
    real(dp) :: minus_negative_volume, minus_negative_centroid(3)
    real(dp) :: minus_negative_second(3, 3), minus_area, minus_normal(3)
    real(dp) :: fd_second_dot(3, 3), fd_centroid_dot(3), fd_volume_dot
    real(dp) :: fd_area_dot, fd_normal_dot(3), step
    integer :: status, plus_status, minus_status

    call evaluate_level_set_tetra_cut_moments_3d( &
        vertices, level_values, positive_volume, positive_centroid, positive_second, &
        negative_volume, negative_centroid, negative_second, interface_area, normal, &
        status)
    call check_condition(status == 0 .and. abs(positive_volume - 1.0_dp/48.0_dp) < &
        1.0e-14_dp .and. maxval(abs(positive_centroid - expected_centroid)) < &
        1.0e-14_dp .and. maxval(abs(positive_second - expected_second)) < &
        1.0e-14_dp, "tetra cut moments match the independent one-corner oracle")
    call check_condition(abs(positive_volume + negative_volume - 1.0_dp/6.0_dp) < &
        1.0e-14_dp .and. maxval(abs(positive_second + negative_second - &
        parent_second)) < 1.0e-14_dp, &
        "tetra cut moments conserve parent volume and quadratic moments")

    call evaluate_level_set_tetra_cut_moments_3d_jvp( &
        vertices, level_values, vertices_dot, level_values_dot, positive_volume_dot, &
        positive_centroid_dot, positive_second_dot, negative_volume_dot, &
        negative_centroid_dot, negative_second_dot, interface_area_dot, normal_dot, &
        status)
    step = 1.0e-6_dp
    call evaluate_level_set_tetra_cut_moments_3d( &
        vertices + step*vertices_dot, level_values + step*level_values_dot, &
        plus_volume, plus_centroid, plus_second, plus_negative_volume, &
        plus_negative_centroid, plus_negative_second, plus_area, plus_normal, &
        plus_status)
    call evaluate_level_set_tetra_cut_moments_3d( &
        vertices - step*vertices_dot, level_values - step*level_values_dot, &
        minus_volume, minus_centroid, minus_second, minus_negative_volume, &
        minus_negative_centroid, minus_negative_second, minus_area, minus_normal, &
        minus_status)
    fd_volume_dot = (plus_volume - minus_volume)/(2.0_dp*step)
    fd_centroid_dot = (plus_centroid - minus_centroid)/(2.0_dp*step)
    fd_second_dot = (plus_second - minus_second)/(2.0_dp*step)
    fd_area_dot = (plus_area - minus_area)/(2.0_dp*step)
    fd_normal_dot = (plus_normal - minus_normal)/(2.0_dp*step)
    call check_condition(status == 0 .and. plus_status == 0 .and. minus_status == 0 .and. &
        abs(positive_volume_dot - fd_volume_dot) < 1.0e-7_dp .and. &
        maxval(abs(positive_centroid_dot - fd_centroid_dot)) < 1.0e-7_dp .and. &
        maxval(abs(positive_second_dot - fd_second_dot)) < 1.0e-7_dp .and. &
        abs(interface_area_dot - fd_area_dot) < 1.0e-7_dp .and. &
        maxval(abs(normal_dot - fd_normal_dot)) < 1.0e-7_dp, &
        "tetra cut-moment JVP matches a fixed-topology central-difference oracle")

    call check_summary("3D level-set tetra cut moments")
end program test_level_set_tetra_cut_moments_3d
