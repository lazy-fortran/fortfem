program test_level_set_triangle_cut_quadrature_jvp_2d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_level_set_triangle_cut_quadrature_2d, &
        evaluate_level_set_triangle_cut_quadrature_2d_jvp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(2, 3) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
    real(dp), parameter :: level_values(3) = [-0.4_dp, 0.6_dp, 0.6_dp]
    real(dp), parameter :: vertices_dot(2, 3) = reshape([ &
        0.03_dp, -0.02_dp, -0.01_dp, 0.04_dp, 0.02_dp, 0.01_dp], [2, 3])
    real(dp), parameter :: level_values_dot(3) = [0.1_dp, -0.03_dp, 0.08_dp]
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: positive_area, positive_area_dot, positive_area_plus
    real(dp) :: positive_area_minus, negative_area, negative_area_dot
    real(dp) :: negative_area_plus, negative_area_minus
    real(dp) :: positive_centroid(2), positive_centroid_dot(2)
    real(dp) :: positive_centroid_plus(2), positive_centroid_minus(2)
    real(dp) :: negative_centroid(2), negative_centroid_dot(2)
    real(dp) :: negative_centroid_plus(2), negative_centroid_minus(2)
    real(dp) :: interface_length, interface_length_dot
    real(dp) :: interface_length_plus, interface_length_minus
    real(dp) :: normal(2), normal_dot(2), normal_plus(2), normal_minus(2)
    real(dp) :: total_area, total_area_dot, total_first_moment(2)
    real(dp) :: total_first_moment_dot(2), first_moment(2), first_moment_dot(2)
    real(dp) :: uncut_area_dot, uncut_centroid_dot(2)
    real(dp) :: uncut_area, uncut_centroid(2)
    integer :: status, status_plus, status_minus

    call evaluate_level_set_triangle_cut_quadrature_2d( &
        vertices, level_values, positive_area, positive_centroid, &
        negative_area, negative_centroid, interface_length, normal, status)
    call check_condition(status == 0, &
        "cut quadrature JVP accepts a proper fixed-topology cut")
    call evaluate_level_set_triangle_cut_quadrature_2d_jvp( &
        vertices, level_values, vertices_dot, level_values_dot, &
        positive_area_dot, positive_centroid_dot, negative_area_dot, &
        negative_centroid_dot, interface_length_dot, normal_dot, status)
    call check_condition(status == 0, &
        "cut quadrature JVP accepts fixed-topology geometry data")

    call evaluate_level_set_triangle_cut_quadrature_2d( &
        vertices + step*vertices_dot, level_values + step*level_values_dot, &
        positive_area_plus, positive_centroid_plus, negative_area_plus, &
        negative_centroid_plus, interface_length_plus, normal_plus, status_plus)
    call evaluate_level_set_triangle_cut_quadrature_2d( &
        vertices - step*vertices_dot, level_values - step*level_values_dot, &
        positive_area_minus, positive_centroid_minus, negative_area_minus, &
        negative_centroid_minus, interface_length_minus, normal_minus, status_minus)
    call check_condition(status_plus == 0 .and. status_minus == 0, &
        "cut quadrature central-difference states retain topology")
    call check_condition(abs(positive_area_dot - &
        (positive_area_plus - positive_area_minus)/(2.0_dp*step)) < 4.0e-9_dp .and. &
        maxval(abs(positive_centroid_dot - &
        (positive_centroid_plus - positive_centroid_minus)/(2.0_dp*step))) < 4.0e-9_dp .and. &
        abs(negative_area_dot - &
        (negative_area_plus - negative_area_minus)/(2.0_dp*step)) < 4.0e-9_dp .and. &
        maxval(abs(negative_centroid_dot - &
        (negative_centroid_plus - negative_centroid_minus)/(2.0_dp*step))) < 4.0e-9_dp .and. &
        abs(interface_length_dot - &
        (interface_length_plus - interface_length_minus)/(2.0_dp*step)) < 4.0e-9_dp .and. &
        maxval(abs(normal_dot - &
        (normal_plus - normal_minus)/(2.0_dp*step))) < 4.0e-9_dp, &
        "cut quadrature JVP matches an independent central difference")

    call triangle_moments( &
        vertices, vertices_dot, total_area, total_area_dot, &
        total_first_moment, total_first_moment_dot)
    first_moment = positive_area*positive_centroid + &
        negative_area*negative_centroid
    first_moment_dot = positive_area_dot*positive_centroid + &
        positive_area*positive_centroid_dot + negative_area_dot*negative_centroid + &
        negative_area*negative_centroid_dot
    call check_condition(abs(positive_area + negative_area - total_area) < 1.0e-14_dp .and. &
        abs(positive_area_dot + negative_area_dot - total_area_dot) < 1.0e-13_dp .and. &
        maxval(abs(first_moment - total_first_moment)) < 1.0e-14_dp .and. &
        maxval(abs(first_moment_dot - total_first_moment_dot)) < 1.0e-13_dp, &
        "cut quadrature JVP preserves area and first-moment conservation")

    call evaluate_level_set_triangle_cut_quadrature_2d_jvp( &
        vertices, [0.2_dp, 0.6_dp, 0.8_dp], vertices_dot, level_values_dot, &
        positive_area_dot, positive_centroid_dot, negative_area_dot, &
        negative_centroid_dot, interface_length_dot, normal_dot, status)
    call triangle_moments(vertices, vertices_dot, uncut_area, uncut_area_dot, &
        total_first_moment, total_first_moment_dot)
    uncut_centroid = sum(vertices, dim=2)/3.0_dp
    uncut_centroid_dot = sum(vertices_dot, dim=2)/3.0_dp
    call check_condition(status == 0 .and. abs(positive_area_dot - uncut_area_dot) < &
        1.0e-13_dp .and. maxval(abs(positive_centroid_dot - uncut_centroid_dot)) < &
        1.0e-13_dp .and. negative_area_dot == 0.0_dp .and. &
        maxval(abs(negative_centroid_dot)) == 0.0_dp .and. &
        interface_length_dot == 0.0_dp .and. maxval(abs(normal_dot)) == 0.0_dp, &
        "cut quadrature JVP handles an uncut positive triangle")
    call evaluate_level_set_triangle_cut_quadrature_2d_jvp( &
        vertices, [0.0_dp, 0.6_dp, 0.8_dp], vertices_dot, level_values_dot, &
        positive_area_dot, positive_centroid_dot, negative_area_dot, &
        negative_centroid_dot, interface_length_dot, normal_dot, status)
    call check_condition(status /= 0, &
        "cut quadrature JVP rejects a nodal topology event")
    call check_summary("2D level-set cut quadrature JVP")

contains

    subroutine triangle_moments( &
            triangle, triangle_dot, area, area_dot, first_moment, &
            first_moment_dot)
        real(dp), intent(in) :: triangle(2, 3), triangle_dot(2, 3)
        real(dp), intent(out) :: area, area_dot, first_moment(2)
        real(dp), intent(out) :: first_moment_dot(2)

        real(dp) :: edge_one(2), edge_two(2), edge_one_dot(2), edge_two_dot(2)
        real(dp) :: signed_area, signed_area_dot, orientation

        edge_one = triangle(:, 2) - triangle(:, 1)
        edge_two = triangle(:, 3) - triangle(:, 1)
        edge_one_dot = triangle_dot(:, 2) - triangle_dot(:, 1)
        edge_two_dot = triangle_dot(:, 3) - triangle_dot(:, 1)
        signed_area = 0.5_dp*(edge_one(1)*edge_two(2) - &
            edge_one(2)*edge_two(1))
        signed_area_dot = 0.5_dp*(edge_one_dot(1)*edge_two(2) + &
            edge_one(1)*edge_two_dot(2) - edge_one_dot(2)*edge_two(1) - &
            edge_one(2)*edge_two_dot(1))
        orientation = sign(1.0_dp, signed_area)
        area = orientation*signed_area
        area_dot = orientation*signed_area_dot
        first_moment = signed_area*sum(triangle, dim=2)/3.0_dp
        first_moment_dot = (signed_area_dot*sum(triangle, dim=2) + &
            signed_area*sum(triangle_dot, dim=2))/3.0_dp
        first_moment = orientation*first_moment
        first_moment_dot = orientation*first_moment_dot
    end subroutine triangle_moments

end program test_level_set_triangle_cut_quadrature_jvp_2d
