program test_level_set_tetra_cut_quadrature_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_level_set_tetra_cut_quadrature_3d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(3, 4) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
    real(dp), parameter :: level_values(4) = [-0.25_dp, 0.75_dp, 0.75_dp, 0.75_dp]
    real(dp), parameter :: total_volume = 1.0_dp/6.0_dp
    real(dp), parameter :: expected_negative_volume = 1.0_dp/384.0_dp
    real(dp), parameter :: expected_negative_centroid(3) = &
        [1.0_dp/16.0_dp, 1.0_dp/16.0_dp, 1.0_dp/16.0_dp]
    real(dp), parameter :: expected_positive_volume = 21.0_dp/128.0_dp
    real(dp), parameter :: expected_positive_centroid(3) = &
        [85.0_dp/336.0_dp, 85.0_dp/336.0_dp, 85.0_dp/336.0_dp]
    real(dp), parameter :: expected_interface_area = sqrt(3.0_dp)/32.0_dp
    real(dp), parameter :: expected_normal(3) = &
        [1.0_dp/sqrt(3.0_dp), 1.0_dp/sqrt(3.0_dp), 1.0_dp/sqrt(3.0_dp)]
    real(dp) :: positive_volume, negative_volume
    real(dp) :: positive_centroid(3), negative_centroid(3)
    real(dp) :: interface_area, normal(3), first_moment(3)
    integer :: status

    call evaluate_level_set_tetra_cut_quadrature_3d( &
        vertices, level_values, positive_volume, positive_centroid, &
        negative_volume, negative_centroid, interface_area, normal, status)
    call check_condition(status == 0, &
        "tetrahedral cut quadrature accepts a proper linear cut")
    call check_condition(abs(positive_volume - expected_positive_volume) < 1.0e-14_dp .and. &
        abs(negative_volume - expected_negative_volume) < 1.0e-14_dp, &
        "tetrahedral cut quadrature matches independent subcell volumes")
    call check_condition(maxval(abs(positive_centroid - expected_positive_centroid)) < &
        1.0e-14_dp .and. maxval(abs(negative_centroid - &
        expected_negative_centroid)) < 1.0e-14_dp, &
        "tetrahedral cut quadrature matches independent subcell centroids")
    call check_condition(abs(interface_area - expected_interface_area) < 1.0e-14_dp .and. &
        maxval(abs(normal - expected_normal)) < 1.0e-14_dp, &
        "tetrahedral cut quadrature preserves interface geometry")
    first_moment = positive_volume*positive_centroid + &
        negative_volume*negative_centroid
    call check_condition(abs(positive_volume + negative_volume - total_volume) < 1.0e-14_dp .and. &
        maxval(abs(first_moment - total_volume/4.0_dp)) < 1.0e-14_dp, &
        "tetrahedral cut quadrature preserves volume and first moment")

    call evaluate_level_set_tetra_cut_quadrature_3d( &
        vertices, [0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp], positive_volume, &
        positive_centroid, negative_volume, negative_centroid, interface_area, &
        normal, status)
    call check_condition(status == 0 .and. abs(positive_volume - total_volume) < &
        1.0e-14_dp .and. maxval(abs(positive_centroid - 0.25_dp)) < &
        1.0e-14_dp .and. negative_volume == 0.0_dp .and. &
        maxval(abs(negative_centroid)) == 0.0_dp .and. interface_area == 0.0_dp .and. &
        maxval(abs(normal)) == 0.0_dp, &
        "tetrahedral cut quadrature handles an uncut positive tetrahedron")

    call evaluate_level_set_tetra_cut_quadrature_3d( &
        vertices, [-0.2_dp, -0.4_dp, -0.6_dp, -0.8_dp], positive_volume, &
        positive_centroid, negative_volume, negative_centroid, interface_area, &
        normal, status)
    call check_condition(status == 0 .and. positive_volume == 0.0_dp .and. &
        maxval(abs(positive_centroid)) == 0.0_dp .and. &
        abs(negative_volume - total_volume) < 1.0e-14_dp .and. &
        maxval(abs(negative_centroid - 0.25_dp)) < 1.0e-14_dp .and. &
        interface_area == 0.0_dp .and. maxval(abs(normal)) == 0.0_dp, &
        "tetrahedral cut quadrature handles an uncut negative tetrahedron")

    call evaluate_level_set_tetra_cut_quadrature_3d( &
        vertices, [-0.5_dp, -0.5_dp, 0.5_dp, 0.5_dp], positive_volume, &
        positive_centroid, negative_volume, negative_centroid, interface_area, &
        normal, status)
    first_moment = positive_volume*positive_centroid + &
        negative_volume*negative_centroid
    call check_condition(status == 0 .and. positive_volume > 0.0_dp .and. &
        negative_volume > 0.0_dp .and. abs(positive_volume + negative_volume - &
        total_volume) < 1.0e-14_dp .and. maxval(abs(first_moment - &
        total_volume/4.0_dp)) < 1.0e-14_dp, &
        "tetrahedral cut quadrature conserves a quadrilateral cut")

    call evaluate_level_set_tetra_cut_quadrature_3d( &
        vertices, [0.0_dp, 0.75_dp, 0.75_dp, 0.75_dp], positive_volume, &
        positive_centroid, negative_volume, negative_centroid, interface_area, &
        normal, status)
    call check_condition(status /= 0, &
        "tetrahedral cut quadrature rejects a nodal topology event")

    call evaluate_level_set_tetra_cut_quadrature_3d( &
        reshape([0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        2.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, 0.0_dp], [3, 4]), &
        level_values, positive_volume, positive_centroid, negative_volume, &
        negative_centroid, interface_area, normal, status)
    call check_condition(status /= 0, &
        "tetrahedral cut quadrature rejects a degenerate map")
    call check_summary("3D level-set tetrahedron cut quadrature")
end program test_level_set_tetra_cut_quadrature_3d
