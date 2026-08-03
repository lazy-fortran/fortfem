program test_level_set_triangle_interface_jvp_2d
    use check, only: check_condition, check_summary
    use fortfem_core, only: &
        evaluate_level_set_triangle_interface_2d, &
        evaluate_level_set_triangle_interface_2d_jvp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(2, 3) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
    real(dp), parameter :: level_values(3) = [-0.4_dp, 0.6_dp, 0.6_dp]
    real(dp), parameter :: vertices_dot(2, 3) = reshape([ &
        0.03_dp, -0.02_dp, -0.01_dp, 0.04_dp, 0.02_dp, 0.01_dp], [2, 3])
    real(dp), parameter :: level_values_dot(3) = [0.1_dp, -0.03_dp, 0.08_dp]
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: points(2, 2), points_dot(2, 2), points_plus(2, 2)
    real(dp) :: points_minus(2, 2), length, length_dot, length_plus
    real(dp) :: length_minus, normal(2), normal_dot(2), normal_plus(2)
    real(dp) :: normal_minus(2), expected_points_dot(2, 2)
    real(dp) :: expected_length_dot, expected_normal_dot(2)
    real(dp) :: edge(2), edge_dot(2), fraction, fraction_dot, denominator
    real(dp) :: gradient(2), gradient_dot(2), gradient_norm, gradient_norm_dot
    integer :: status, status_plus, status_minus

    call evaluate_level_set_triangle_interface_2d( &
        vertices, level_values, points, length, normal, status)
    call check_condition(status == 0, &
        "level-set interface JVP accepts a fixed proper cut")
    call evaluate_level_set_triangle_interface_2d_jvp( &
        vertices, level_values, vertices_dot, level_values_dot, points_dot, &
        length_dot, normal_dot, status)
    call check_condition(status == 0, &
        "level-set interface JVP accepts fixed-topology geometry data")

    call intersection_dot(1, 2, expected_points_dot(:, 1))
    call intersection_dot(1, 3, expected_points_dot(:, 2))
    call expected_gradient_and_dot( &
        gradient, gradient_dot, gradient_norm, gradient_norm_dot)
    edge = points(:, 2) - points(:, 1)
    edge_dot = expected_points_dot(:, 2) - expected_points_dot(:, 1)
    expected_length_dot = dot_product(edge, edge_dot)/length
    expected_normal_dot = (gradient_dot*gradient_norm - &
        gradient*gradient_norm_dot)/(gradient_norm*gradient_norm)
    call check_condition(maxval(abs(points_dot - expected_points_dot)) < 1.0e-13_dp .and. &
        abs(length_dot - expected_length_dot) < 1.0e-13_dp .and. &
        maxval(abs(normal_dot - expected_normal_dot)) < 1.0e-13_dp, &
        "level-set interface JVP matches independent intersection and normal formulas")

    call evaluate_level_set_triangle_interface_2d( &
        vertices + step*vertices_dot, level_values + step*level_values_dot, &
        points_plus, length_plus, normal_plus, status_plus)
    call evaluate_level_set_triangle_interface_2d( &
        vertices - step*vertices_dot, level_values - step*level_values_dot, &
        points_minus, length_minus, normal_minus, status_minus)
    call check_condition(status_plus == 0 .and. status_minus == 0, &
        "level-set interface central-difference states retain topology")
    call check_condition(maxval(abs(points_dot - &
        (points_plus - points_minus)/(2.0_dp*step))) < 3.0e-9_dp .and. &
        abs(length_dot - (length_plus - length_minus)/(2.0_dp*step)) < 3.0e-9_dp .and. &
        maxval(abs(normal_dot - &
        (normal_plus - normal_minus)/(2.0_dp*step))) < 3.0e-9_dp, &
        "level-set interface JVP matches central differences")

    call evaluate_level_set_triangle_interface_2d_jvp( &
        vertices, [0.0_dp, 0.6_dp, 0.6_dp], vertices_dot, level_values_dot, &
        points_dot, length_dot, normal_dot, status)
    call check_condition(status /= 0, &
        "level-set interface JVP rejects a vertex topology event")
    call check_summary("2D level-set triangle interface JVP")

contains

    subroutine intersection_dot(first, second, point_dot)
        integer, intent(in) :: first, second
        real(dp), intent(out) :: point_dot(2)

        real(dp) :: first_value, second_value, first_dot, second_dot
        real(dp) :: edge_local(2), edge_local_dot(2), local_fraction
        real(dp) :: local_fraction_dot, local_denominator

        first_value = level_values(first)
        second_value = level_values(second)
        first_dot = level_values_dot(first)
        second_dot = level_values_dot(second)
        local_denominator = first_value - second_value
        local_fraction = first_value/local_denominator
        local_fraction_dot = (first_dot*local_denominator - &
            first_value*(first_dot - second_dot))/(local_denominator**2)
        edge_local = vertices(:, second) - vertices(:, first)
        edge_local_dot = vertices_dot(:, second) - vertices_dot(:, first)
        point_dot = vertices_dot(:, first) + local_fraction_dot*edge_local + &
            local_fraction*edge_local_dot
    end subroutine intersection_dot

    subroutine expected_gradient_and_dot( &
            expected_gradient, expected_gradient_dot, norm, norm_dot)
        real(dp), intent(out) :: expected_gradient(2), expected_gradient_dot(2)
        real(dp), intent(out) :: norm, norm_dot

        real(dp) :: first_edge(2), second_edge(2), first_edge_dot(2)
        real(dp) :: second_edge_dot(2), determinant, determinant_dot
        real(dp) :: level_a, level_b, level_a_dot, level_b_dot
        real(dp) :: numerator_x, numerator_y, numerator_x_dot, numerator_y_dot

        first_edge = vertices(:, 2) - vertices(:, 1)
        second_edge = vertices(:, 3) - vertices(:, 1)
        first_edge_dot = vertices_dot(:, 2) - vertices_dot(:, 1)
        second_edge_dot = vertices_dot(:, 3) - vertices_dot(:, 1)
        determinant = first_edge(1)*second_edge(2) - &
            first_edge(2)*second_edge(1)
        determinant_dot = first_edge_dot(1)*second_edge(2) + &
            first_edge(1)*second_edge_dot(2) - first_edge_dot(2)*second_edge(1) - &
            first_edge(2)*second_edge_dot(1)
        level_a = level_values(2) - level_values(1)
        level_b = level_values(3) - level_values(1)
        level_a_dot = level_values_dot(2) - level_values_dot(1)
        level_b_dot = level_values_dot(3) - level_values_dot(1)
        numerator_x = level_a*second_edge(2) - level_b*first_edge(2)
        numerator_y = first_edge(1)*level_b - second_edge(1)*level_a
        numerator_x_dot = level_a_dot*second_edge(2) + &
            level_a*second_edge_dot(2) - level_b_dot*first_edge(2) - &
            level_b*first_edge_dot(2)
        numerator_y_dot = first_edge_dot(1)*level_b + &
            first_edge(1)*level_b_dot - second_edge_dot(1)*level_a - &
            second_edge(1)*level_a_dot
        expected_gradient = [numerator_x, numerator_y]/determinant
        expected_gradient_dot = [ &
            (numerator_x_dot*determinant - numerator_x*determinant_dot)/ &
            determinant**2, &
            (numerator_y_dot*determinant - numerator_y*determinant_dot)/ &
            determinant**2]
        norm = sqrt(dot_product(expected_gradient, expected_gradient))
        norm_dot = dot_product(expected_gradient, expected_gradient_dot)/norm
    end subroutine expected_gradient_and_dot

end program test_level_set_triangle_interface_jvp_2d
