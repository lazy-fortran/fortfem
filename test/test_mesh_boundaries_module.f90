program test_mesh_boundaries_module
    use fortfem_kinds, only: dp
    use fortfem_boundary, only: boundary_t
    use fortfem_api_mesh_boundaries, only: circle_boundary,                &
        rectangle_boundary, line_segment, arc_segment, l_shape_boundary
    use check, only: check_condition, check_summary
    implicit none

    type(boundary_t) :: circle, square, line, arc, lshape
    real(dp) :: radius, radius_point
    integer :: n, n_label_1, n_label_2, n_label_3, n_label_4

    write(*,*) "Testing fortfem_api_mesh_boundaries module..."

    ! Circle boundary: basic geometry and closure
    radius = 1.0_dp
    n = 16
    circle = circle_boundary([0.0_dp, 0.0_dp], radius, n)

    call check_condition(circle%n_points == n,                             &
        "circle_boundary: correct number of points")
    call check_condition(size(circle%labels) == n - 1,                     &
        "circle_boundary: correct number of edge labels")
    call check_condition(circle%is_closed,                                  &
        "circle_boundary: boundary is closed")

    radius_point = sqrt(circle%points(1, 1)**2 + circle%points(2, 1)**2)
    call check_condition(abs(radius_point - radius) < 1.0e-10_dp,          &
        "circle_boundary: points lie on circle")

    ! Rectangle boundary: label distribution and closure
    n = 5
    square = rectangle_boundary([0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], n)

    call check_condition(square%n_points == 4 * n,                          &
        "rectangle_boundary: correct number of points")
    call check_condition(size(square%labels) == 4 * n - 1,                  &
        "rectangle_boundary: correct number of edge labels")
    call check_condition(square%is_closed,                                   &
        "rectangle_boundary: boundary is closed")

    n_label_1 = count(square%labels == 1)
    n_label_2 = count(square%labels == 2)
    n_label_3 = count(square%labels == 3)
    n_label_4 = count(square%labels == 4)

    call check_condition(n_label_1 == n - 1,                                &
        "rectangle_boundary: correct count for label 1")
    call check_condition(n_label_2 == n - 1,                                &
        "rectangle_boundary: correct count for label 2")
    call check_condition(n_label_3 == n - 1,                                &
        "rectangle_boundary: correct count for label 3")
    call check_condition(n_label_4 == n + 2,                                &
        "rectangle_boundary: correct count for label 4")

    ! Line segment: open boundary with uniform spacing
    line = line_segment([0.0_dp, 0.0_dp], [1.0_dp, 0.0_dp], 4)

    call check_condition(.not. line%is_closed,                              &
        "line_segment: boundary is open")
    call check_condition(line%n_points == 4,                                &
        "line_segment: correct number of points")

    ! Arc segment: open boundary between two points on circle
    arc = arc_segment([1.0_dp, 0.0_dp], [0.0_dp, 1.0_dp],                  &
        [0.0_dp, 0.0_dp], 5)

    call check_condition(.not. arc%is_closed,                               &
        "arc_segment: boundary is open")
    call check_condition(arc%n_points == 5,                                 &
        "arc_segment: correct number of points")

    ! L-shape boundary: closed boundary with expected point count
    lshape = l_shape_boundary(1.0_dp, 24)

    call check_condition(lshape%is_closed,                                  &
        "l_shape_boundary: boundary is closed")
    call check_condition(lshape%n_points > 0,                               &
        "l_shape_boundary: positive number of points")

    call check_summary("fortfem_api_mesh_boundaries")

end program test_mesh_boundaries_module

