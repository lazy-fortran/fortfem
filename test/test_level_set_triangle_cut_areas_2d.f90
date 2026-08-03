program test_level_set_triangle_cut_areas_2d
    use check, only: check_condition, check_summary
    use fortfem_core, only: evaluate_level_set_triangle_cut_areas_2d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(2, 3) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
    real(dp), parameter :: level_values(3) = [-0.4_dp, 0.6_dp, 0.6_dp]
    real(dp) :: positive_area, negative_area, interface_length
    integer :: status

    call evaluate_level_set_triangle_cut_areas_2d( &
        vertices, level_values, positive_area, negative_area, &
        interface_length, status)
    call check_condition(status == 0, &
        "level-set cut areas accept a proper linear cut")
    call check_condition(abs(positive_area - 0.42_dp) < 1.0e-14_dp .and. &
        abs(negative_area - 0.08_dp) < 1.0e-14_dp, &
        "level-set cut areas match the independent clipped-triangle oracle")
    call check_condition(abs(interface_length - sqrt(0.32_dp)) < 1.0e-14_dp, &
        "level-set cut area interface length matches the segment oracle")
    call check_condition(abs(positive_area + negative_area - 0.5_dp) < &
        1.0e-14_dp, "level-set cut areas conserve the parent triangle area")

    call evaluate_level_set_triangle_cut_areas_2d( &
        vertices, [0.2_dp, 0.6_dp, 0.8_dp], positive_area, negative_area, &
        interface_length, status)
    call check_condition(status == 0 .and. abs(positive_area - 0.5_dp) < &
        1.0e-14_dp .and. negative_area == 0.0_dp .and. &
        interface_length == 0.0_dp, &
        "level-set cut areas handle an uncut positive triangle")
    call check_summary("2D level-set triangle cut areas")
end program test_level_set_triangle_cut_areas_2d
