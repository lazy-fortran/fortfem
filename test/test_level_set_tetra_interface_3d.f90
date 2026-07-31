program test_level_set_tetra_interface_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_level_set_tetra_interface_3d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(3, 4) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
    real(dp), parameter :: level_values(4) = [-0.25_dp, 0.75_dp, 0.75_dp, 0.75_dp]
    real(dp), parameter :: expected_points(3, 3) = reshape([ &
        0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.25_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.25_dp], [3, 3])
    real(dp), parameter :: expected_normal(3) = &
        [1.0_dp/sqrt(3.0_dp), 1.0_dp/sqrt(3.0_dp), 1.0_dp/sqrt(3.0_dp)]
    real(dp), parameter :: expected_area = sqrt(3.0_dp)/32.0_dp
    real(dp) :: points(3, 4), area, normal(3)
    real(dp) :: distances(3), expected_area_vector(3)
    integer :: point_count, status

    call evaluate_level_set_tetra_interface_3d( &
        vertices, level_values, points, point_count, area, normal, status)
    call check_condition(status == 0 .and. point_count == 3, &
        "tetrahedral level-set interface accepts a triangular cut")
    call check_condition(maxval(abs(normal - expected_normal)) < 1.0e-14_dp .and. &
        abs(area - expected_area) < 1.0e-14_dp, &
        "tetrahedral interface matches the independent plane normal and area")
    call match_points(points(:, :point_count), expected_points, distances)
    call check_condition(maxval(distances) < 1.0e-14_dp, &
        "tetrahedral interface returns all independent edge intersections")

    expected_area_vector = [1.0_dp, 1.0_dp, 1.0_dp]/32.0_dp
    call check_condition(abs(dot_product(normal, expected_area_vector) - &
        3.0_dp/(32.0_dp*sqrt(3.0_dp))) < 1.0e-14_dp, &
        "tetrahedral interface preserves the oriented gradient convention")

    call evaluate_level_set_tetra_interface_3d( &
        vertices, [-0.5_dp, -0.5_dp, 0.5_dp, 0.5_dp], points, point_count, &
        area, normal, status)
    call check_condition(status == 0 .and. point_count == 4 .and. area > 0.0_dp, &
        "tetrahedral level-set interface accepts a quadrilateral cut")

    call evaluate_level_set_tetra_interface_3d( &
        vertices, [0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp], points, point_count, &
        area, normal, status)
    call check_condition(status /= 0 .and. point_count == 0, &
        "tetrahedral level-set interface rejects an uncut tetrahedron")

    call evaluate_level_set_tetra_interface_3d( &
        vertices, [0.0_dp, 0.75_dp, 0.75_dp, 0.75_dp], points, point_count, &
        area, normal, status)
    call check_condition(status /= 0, &
        "tetrahedral level-set interface rejects a nodal topology event")

    call evaluate_level_set_tetra_interface_3d( &
        reshape([0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        2.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, 0.0_dp], [3, 4]), &
        level_values, points, point_count, area, normal, status)
    call check_condition(status /= 0, &
        "tetrahedral level-set interface rejects a degenerate map")
    call check_summary("3D level-set tetrahedron interface")

contains

    subroutine match_points(actual, expected, distances)
        real(dp), intent(in) :: actual(:, :), expected(:, :)
        real(dp), intent(out) :: distances(:)

        integer :: point, candidate

        distances = huge(1.0_dp)
        do point = 1, size(expected, 2)
            do candidate = 1, size(actual, 2)
                distances(point) = min(distances(point), &
                    sqrt(sum((expected(:, point) - actual(:, candidate))**2)))
            end do
        end do
    end subroutine match_points

end program test_level_set_tetra_interface_3d
