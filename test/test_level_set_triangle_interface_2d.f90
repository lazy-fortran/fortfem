program test_level_set_triangle_interface_2d
    use check, only: check_condition, check_summary
    use fortfem_core, only: evaluate_level_set_triangle_interface_2d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(2, 3) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
    real(dp), parameter :: level_values(3) = [-0.4_dp, 0.6_dp, 0.6_dp]
    real(dp), parameter :: expected_points(2, 2) = reshape([ &
        0.4_dp, 0.0_dp, 0.0_dp, 0.4_dp], [2, 2])
    real(dp) :: points(2, 2), normal(2), length
    integer :: status

    call evaluate_level_set_triangle_interface_2d( &
        vertices, level_values, points, length, normal, status)
    call check_condition(status == 0, &
        "level-set triangle interface accepts a proper cut")
    call check_condition(maxval(abs(points - expected_points)) < 1.0e-14_dp, &
        "level-set interface returns the independent edge intersections")
    call check_condition(abs(length - sqrt(0.32_dp)) < 1.0e-14_dp, &
        "level-set interface returns the physical segment length")
    call check_condition(maxval(abs(normal - &
        [1.0_dp/sqrt(2.0_dp), 1.0_dp/sqrt(2.0_dp)])) < 1.0e-14_dp, &
        "level-set interface normal follows the level-set gradient")

    call evaluate_level_set_triangle_interface_2d( &
        vertices, [0.2_dp, 0.6_dp, 0.8_dp], points, length, normal, status)
    call check_condition(status /= 0, &
        "level-set interface rejects an uncut triangle")
    call check_summary("2D level-set triangle interface")
end program test_level_set_triangle_interface_2d
