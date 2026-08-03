program test_level_set_triangle_cut_quadrature_2d
    use check, only: check_condition, check_summary
    use fortfem_core, only: evaluate_level_set_triangle_cut_quadrature_2d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(2, 3) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
    real(dp), parameter :: level_values(3) = [-0.4_dp, 0.6_dp, 0.6_dp]
    real(dp), parameter :: expected_positive_centroid(2) = [ &
        13.0_dp/35.0_dp, 13.0_dp/35.0_dp]
    real(dp), parameter :: expected_negative_centroid(2) = [ &
        2.0_dp/15.0_dp, 2.0_dp/15.0_dp]
    real(dp), parameter :: expected_normal(2) = [ &
        1.0_dp/sqrt(2.0_dp), 1.0_dp/sqrt(2.0_dp)]
    real(dp) :: positive_area, negative_area
    real(dp) :: positive_centroid(2), negative_centroid(2)
    real(dp) :: interface_length, normal(2), positive_affine_integral
    integer :: status

    call evaluate_level_set_triangle_cut_quadrature_2d( &
        vertices, level_values, positive_area, positive_centroid, &
        negative_area, negative_centroid, interface_length, normal, status)
    call check_condition(status == 0, &
        "level-set cut quadrature accepts a proper linear cut")
    call check_condition(abs(positive_area - 0.42_dp) < 1.0e-14_dp .and. &
        abs(negative_area - 0.08_dp) < 1.0e-14_dp, &
        "level-set cut quadrature preserves exact subcell areas")
    call check_condition(maxval(abs(positive_centroid - &
        expected_positive_centroid)) < 1.0e-14_dp .and. &
        maxval(abs(negative_centroid - expected_negative_centroid)) < &
        1.0e-14_dp, "level-set cut quadrature matches polygon centroids")
    call check_condition(abs(interface_length - sqrt(0.32_dp)) < 1.0e-14_dp .and. &
        maxval(abs(normal - expected_normal)) < 1.0e-14_dp, &
        "level-set cut quadrature preserves interface geometry")
    positive_affine_integral = positive_area*( &
        2.0_dp*positive_centroid(1) - positive_centroid(2) + 3.0_dp)
    call check_condition(abs(positive_affine_integral - 1.416_dp) < 1.0e-14_dp, &
        "level-set cut quadrature integrates affine fields exactly")

    call evaluate_level_set_triangle_cut_quadrature_2d( &
        vertices, [0.2_dp, 0.6_dp, 0.8_dp], positive_area, positive_centroid, &
        negative_area, negative_centroid, interface_length, normal, status)
    call check_condition(status == 0 .and. abs(positive_area - 0.5_dp) < &
        1.0e-14_dp .and. maxval(abs(positive_centroid - &
        [1.0_dp/3.0_dp, 1.0_dp/3.0_dp])) < 1.0e-14_dp .and. &
        negative_area == 0.0_dp .and. maxval(abs(negative_centroid)) == 0.0_dp .and. &
        interface_length == 0.0_dp .and. maxval(abs(normal)) == 0.0_dp, &
        "level-set cut quadrature handles an uncut positive triangle")

    call evaluate_level_set_triangle_cut_quadrature_2d( &
        reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp], [2, 3]), &
        level_values, positive_area, positive_centroid, negative_area, &
        negative_centroid, interface_length, normal, status)
    call check_condition(status /= 0, &
        "level-set cut quadrature rejects a degenerate parent triangle")
    call check_summary("2D level-set triangle cut quadrature")
end program test_level_set_triangle_cut_quadrature_2d
