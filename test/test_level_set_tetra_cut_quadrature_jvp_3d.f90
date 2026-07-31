program test_level_set_tetra_cut_quadrature_jvp_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_level_set_tetra_cut_quadrature_3d, &
        evaluate_level_set_tetra_cut_quadrature_3d_jvp, &
        evaluate_level_set_tetra_interface_3d_jvp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(3, 4) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
    real(dp), parameter :: level_values(4) = [-0.25_dp, 0.75_dp, 0.75_dp, 0.75_dp]
    real(dp), parameter :: vertices_dot(3, 4) = reshape([ &
        0.03_dp, -0.02_dp, 0.01_dp, -0.01_dp, 0.04_dp, 0.02_dp, &
        0.02_dp, -0.03_dp, 0.05_dp, -0.04_dp, 0.01_dp, -0.02_dp], [3, 4])
    real(dp), parameter :: level_values_dot(4) = [0.1_dp, -0.03_dp, 0.08_dp, 0.02_dp]
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: positive_volume, negative_volume
    real(dp) :: positive_centroid(3), negative_centroid(3)
    real(dp) :: interface_area, normal(3)
    real(dp) :: positive_volume_dot, negative_volume_dot
    real(dp) :: positive_centroid_dot(3), negative_centroid_dot(3)
    real(dp) :: interface_area_dot, normal_dot(3)
    real(dp) :: positive_volume_plus, negative_volume_plus
    real(dp) :: positive_centroid_plus(3), negative_centroid_plus(3)
    real(dp) :: interface_area_plus, normal_plus(3)
    real(dp) :: positive_volume_minus, negative_volume_minus
    real(dp) :: positive_centroid_minus(3), negative_centroid_minus(3)
    real(dp) :: interface_area_minus, normal_minus(3)
    real(dp) :: points_dot(3, 4), interface_area_dot_only, normal_dot_only(3)
    real(dp) :: total_volume, total_volume_dot, total_moment(3)
    real(dp) :: total_moment_dot(3), first_moment(3), first_moment_dot(3)
    integer :: status, status_plus, status_minus, interface_status

    call evaluate_level_set_tetra_cut_quadrature_3d( &
        vertices, level_values, positive_volume, positive_centroid, &
        negative_volume, negative_centroid, interface_area, normal, status)
    call check_condition(status == 0, &
        "tetrahedral cut quadrature JVP accepts a proper cut")
    call evaluate_level_set_tetra_cut_quadrature_3d_jvp( &
        vertices, level_values, vertices_dot, level_values_dot, &
        positive_volume_dot, positive_centroid_dot, negative_volume_dot, &
        negative_centroid_dot, interface_area_dot, normal_dot, status)
    call check_condition(status == 0, &
        "tetrahedral cut quadrature JVP accepts fixed topology")

    call evaluate_level_set_tetra_cut_quadrature_3d( &
        vertices + step*vertices_dot, level_values + step*level_values_dot, &
        positive_volume_plus, positive_centroid_plus, negative_volume_plus, &
        negative_centroid_plus, interface_area_plus, normal_plus, status_plus)
    call evaluate_level_set_tetra_cut_quadrature_3d( &
        vertices - step*vertices_dot, level_values - step*level_values_dot, &
        positive_volume_minus, positive_centroid_minus, negative_volume_minus, &
        negative_centroid_minus, interface_area_minus, normal_minus, status_minus)
    call check_condition(status_plus == 0 .and. status_minus == 0, &
        "tetrahedral cut central-difference states retain topology")
    call check_condition(abs(positive_volume_dot - &
        (positive_volume_plus - positive_volume_minus)/(2.0_dp*step)) < 4.0e-9_dp .and. &
        maxval(abs(positive_centroid_dot - (positive_centroid_plus - &
        positive_centroid_minus)/(2.0_dp*step))) < 4.0e-9_dp .and. &
        abs(negative_volume_dot - (negative_volume_plus - &
        negative_volume_minus)/(2.0_dp*step)) < 4.0e-9_dp .and. &
        maxval(abs(negative_centroid_dot - (negative_centroid_plus - &
        negative_centroid_minus)/(2.0_dp*step))) < 4.0e-9_dp .and. &
        abs(interface_area_dot - (interface_area_plus - &
        interface_area_minus)/(2.0_dp*step)) < 4.0e-9_dp .and. &
        maxval(abs(normal_dot - (normal_plus - normal_minus)/(2.0_dp*step))) < &
        4.0e-9_dp, "tetrahedral cut JVP matches an independent central difference")

    call tetrahedron_moments( &
        vertices, vertices_dot, total_volume, total_volume_dot, total_moment, &
        total_moment_dot)
    first_moment = positive_volume*positive_centroid + &
        negative_volume*negative_centroid
    first_moment_dot = positive_volume_dot*positive_centroid + &
        positive_volume*positive_centroid_dot + negative_volume_dot*negative_centroid + &
        negative_volume*negative_centroid_dot
    call check_condition(abs(positive_volume_dot + negative_volume_dot - &
        total_volume_dot) < 1.0e-12_dp .and. maxval(abs(first_moment - &
        total_moment)) < 1.0e-14_dp .and. maxval(abs(first_moment_dot - &
        total_moment_dot)) < 1.0e-12_dp, &
        "tetrahedral cut JVP preserves volume and first-moment conservation")

    call evaluate_level_set_tetra_interface_3d_jvp( &
        vertices, level_values, vertices_dot, level_values_dot, points_dot, &
        interface_area_dot_only, normal_dot_only, interface_status)
    call check_condition(interface_status == 0 .and. &
        abs(interface_area_dot_only - interface_area_dot) < 1.0e-14_dp .and. &
        maxval(abs(normal_dot_only - normal_dot)) < 1.0e-14_dp, &
        "tetrahedral interface JVP shares the cut quadrature contract")

    call evaluate_level_set_tetra_cut_quadrature_3d_jvp( &
        vertices, [0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp], vertices_dot, &
        level_values_dot, positive_volume_dot, positive_centroid_dot, &
        negative_volume_dot, negative_centroid_dot, interface_area_dot, &
        normal_dot, status)
    call check_condition(status == 0 .and. abs(positive_volume_dot - &
        total_volume_dot) < 1.0e-12_dp .and. maxval(abs(positive_centroid_dot - &
        (sum(vertices_dot, dim=2)/4.0_dp))) < 1.0e-12_dp .and. &
        negative_volume_dot == 0.0_dp .and. maxval(abs(negative_centroid_dot)) == 0.0_dp .and. &
        interface_area_dot == 0.0_dp .and. maxval(abs(normal_dot)) == 0.0_dp, &
        "tetrahedral cut JVP handles an uncut positive tetrahedron")

    call evaluate_level_set_tetra_cut_quadrature_3d_jvp( &
        vertices, [0.0_dp, 0.75_dp, 0.75_dp, 0.75_dp], vertices_dot, &
        level_values_dot, positive_volume_dot, positive_centroid_dot, &
        negative_volume_dot, negative_centroid_dot, interface_area_dot, &
        normal_dot, status)
    call check_condition(status /= 0, &
        "tetrahedral cut JVP rejects a nodal topology event")
    call check_summary("3D level-set tetrahedron cut quadrature JVP")

contains

    subroutine tetrahedron_moments( &
            tetrahedron, tetrahedron_dot, volume, volume_dot, moment, moment_dot)
        real(dp), intent(in) :: tetrahedron(3, 4), tetrahedron_dot(3, 4)
        real(dp), intent(out) :: volume, volume_dot, moment(3), moment_dot(3)

        real(dp) :: edge_one(3), edge_two(3), edge_three(3)
        real(dp) :: edge_one_dot(3), edge_two_dot(3), edge_three_dot(3)
        real(dp) :: determinant, determinant_dot

        edge_one = tetrahedron(:, 2) - tetrahedron(:, 1)
        edge_two = tetrahedron(:, 3) - tetrahedron(:, 1)
        edge_three = tetrahedron(:, 4) - tetrahedron(:, 1)
        edge_one_dot = tetrahedron_dot(:, 2) - tetrahedron_dot(:, 1)
        edge_two_dot = tetrahedron_dot(:, 3) - tetrahedron_dot(:, 1)
        edge_three_dot = tetrahedron_dot(:, 4) - tetrahedron_dot(:, 1)
        determinant = dot_product(edge_one, cross_product(edge_two, edge_three))
        determinant_dot = dot_product(edge_one_dot, cross_product(edge_two, edge_three)) + &
            dot_product(edge_one, cross_product(edge_two_dot, edge_three)) + &
            dot_product(edge_one, cross_product(edge_two, edge_three_dot))
        volume = abs(determinant)/6.0_dp
        volume_dot = sign(1.0_dp, determinant)*determinant_dot/6.0_dp
        moment = volume*sum(tetrahedron, dim=2)/4.0_dp
        moment_dot = volume_dot*sum(tetrahedron, dim=2)/4.0_dp + &
            volume*sum(tetrahedron_dot, dim=2)/4.0_dp
    end subroutine tetrahedron_moments

    pure function cross_product(first, second) result(cross)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: cross(3)

        cross = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end program test_level_set_tetra_cut_quadrature_jvp_3d
