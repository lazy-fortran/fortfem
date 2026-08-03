program test_fci_curved_polygon_area
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        compute_fci_curved_polygon_cell_areas_2d, &
        compute_fci_curved_polygon_cell_areas_2d_jvp, &
        compute_fci_curved_polygon_cell_areas_2d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: vertex_count = 5
    integer, parameter :: cell_count = 2
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp), parameter :: gauss_points(3) = [ &
        0.1127016653792583_dp, 0.5_dp, 0.8872983346207417_dp]
    real(dp), parameter :: gauss_weights(3) = [ &
        5.0_dp/18.0_dp, 8.0_dp/18.0_dp, 5.0_dp/18.0_dp]
    real(dp) :: vertices(2, vertex_count, cell_count)
    real(dp) :: controls(2, vertex_count, cell_count)
    real(dp) :: vertices_dot(2, vertex_count, cell_count)
    real(dp) :: controls_dot(2, vertex_count, cell_count)
    real(dp) :: vertices_bar(2, vertex_count, cell_count)
    real(dp) :: controls_bar(2, vertex_count, cell_count)
    real(dp) :: vertices_plus(2, vertex_count, cell_count)
    real(dp) :: vertices_minus(2, vertex_count, cell_count)
    real(dp) :: controls_plus(2, vertex_count, cell_count)
    real(dp) :: controls_minus(2, vertex_count, cell_count)
    real(dp) :: midpoint_controls(2, vertex_count, cell_count)
    real(dp) :: areas(cell_count), areas_dot(cell_count)
    real(dp) :: areas_bar(cell_count), areas_plus(cell_count)
    real(dp) :: areas_minus(cell_count), expected(cell_count)
    real(dp) :: straight(cell_count), polygon_expected(cell_count)
    real(dp) :: lhs, rhs
    real(dp) :: invalid_vertices(2, vertex_count, cell_count)
    type(fortsparse_status_t) :: status
    integer :: cell, edge, next_edge
    logical :: all_passed

    all_passed = .true.
    vertices(:, :, 1) = reshape([ &
        0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 1.0_dp, &
        1.6_dp, 2.0_dp, 0.0_dp, 1.4_dp], [2, vertex_count])
    vertices(:, :, 2) = reshape([ &
        4.0_dp, 0.0_dp, 6.0_dp, 0.3_dp, 6.4_dp, 1.7_dp, &
        5.0_dp, 2.3_dp, 3.9_dp, 1.1_dp], [2, vertex_count])
    controls(:, :, 1) = reshape([ &
        1.0_dp, -0.25_dp, 2.7_dp, 0.5_dp, 2.4_dp, 1.7_dp, &
        0.7_dp, 2.1_dp, -0.2_dp, 0.6_dp], [2, vertex_count])
    controls(:, :, 2) = reshape([ &
        5.0_dp, 0.0_dp, 6.6_dp, 1.0_dp, 5.6_dp, 2.4_dp, &
        4.2_dp, 1.8_dp, 3.8_dp, 0.4_dp], [2, vertex_count])
    do cell = 1, cell_count
        do edge = 1, vertex_count
            next_edge = mod(edge, vertex_count) + 1
            midpoint_controls(:, edge, cell) = 0.5_dp*( &
                vertices(:, edge, cell) + vertices(:, next_edge, cell))
            vertices_dot(:, edge, cell) = [ &
                0.01_dp*sin(real(3*edge + cell, dp)), &
                -0.02_dp*cos(real(2*edge + cell, dp))]
            controls_dot(:, edge, cell) = [ &
                -0.015_dp*cos(real(edge + 2*cell, dp)), &
                0.012_dp*sin(real(2*edge + cell, dp))]
        end do
    end do
    vertices_bar = reshape([(0.03_dp*sin(real(5*edge, dp)), &
        edge=1, size(vertices_bar))], shape(vertices_bar))
    controls_bar = reshape([(0.02_dp*cos(real(3*edge, dp)), &
        edge=1, size(controls_bar))], shape(controls_bar))
    areas_bar = [0.7_dp, -1.1_dp]

    call compute_fci_curved_polygon_cell_areas_2d( &
        vertices, controls, areas, status)
    call record_condition(status%code == 0, &
        "FCI curved polygon area accepts curved pentagons")
    do cell = 1, cell_count
        expected(cell) = gauss_green_area( &
            vertices(:, :, cell), controls(:, :, cell))
        straight(cell) = shoelace_area(vertices(:, :, cell))
        polygon_expected(cell) = straight(cell)
    end do
    call record_condition(maxval(abs(areas - expected)) < 2.0e-12_dp, &
        "FCI curved polygon area matches Gauss Green oracle")
    call record_condition(maxval(abs(areas - straight)) > 1.0e-3_dp, &
        "FCI curved polygon area reflects edge curvature")
    call record_condition(all(areas > 0.0_dp), &
        "FCI curved polygon areas are positive")

    call compute_fci_curved_polygon_cell_areas_2d_jvp( &
        vertices, controls, vertices_dot, controls_dot, areas_dot, status)
    vertices_plus = vertices + step*vertices_dot
    vertices_minus = vertices - step*vertices_dot
    controls_plus = controls + step*controls_dot
    controls_minus = controls - step*controls_dot
    call compute_fci_curved_polygon_cell_areas_2d( &
        vertices_plus, controls_plus, areas_plus, status)
    call compute_fci_curved_polygon_cell_areas_2d( &
        vertices_minus, controls_minus, areas_minus, status)
    call record_condition(maxval(abs(areas_dot - &
        (areas_plus - areas_minus)/(2.0_dp*step))) < 4.0e-8_dp, &
        "FCI curved polygon area JVP matches central differences")

    call compute_fci_curved_polygon_cell_areas_2d_vjp( &
        vertices, controls, areas_bar, vertices_bar, controls_bar, status)
    lhs = dot_product(areas_bar, areas_dot)
    rhs = sum(vertices_bar*vertices_dot) + sum(controls_bar*controls_dot)
    call record_condition(abs(lhs - rhs) < &
        5.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "FCI curved polygon area VJP satisfies the real dot-product identity")

    call compute_fci_curved_polygon_cell_areas_2d( &
        vertices, midpoint_controls, straight, status)
    call record_condition(maxval(abs(straight - polygon_expected)) < &
        2.0e-12_dp, "FCI curved polygon area recovers straight polygon area")

    invalid_vertices = vertices
    invalid_vertices(:, vertex_count, 1) = invalid_vertices(:, 1, 1)
    call compute_fci_curved_polygon_cell_areas_2d( &
        invalid_vertices, controls, areas, status)
    call record_condition(status%code /= 0, &
        "FCI curved polygon area rejects repeated vertices")

    call check_summary("FCI curved polygon area")
    if (.not. all_passed) error stop 1

contains

    pure real(dp) function shoelace_area(cell_vertices) result(area)
        real(dp), intent(in) :: cell_vertices(2, vertex_count)
        integer :: edge, next_edge

        area = 0.0_dp
        do edge = 1, vertex_count
            next_edge = mod(edge, vertex_count) + 1
            area = area + cell_vertices(1, edge)* &
                cell_vertices(2, next_edge) - &
                cell_vertices(2, edge)*cell_vertices(1, next_edge)
        end do
        area = 0.5_dp*area
    end function shoelace_area

    pure real(dp) function gauss_green_area(cell_vertices, edge_controls) &
            result(area)
        real(dp), intent(in) :: cell_vertices(2, vertex_count)
        real(dp), intent(in) :: edge_controls(2, vertex_count)
        real(dp) :: point(2), tangent(2), contribution
        integer :: edge, quadrature, next_edge

        area = 0.0_dp
        do edge = 1, vertex_count
            next_edge = mod(edge, vertex_count) + 1
            do quadrature = 1, 3
                call bezier_point_and_tangent( &
                    cell_vertices(:, edge), edge_controls(:, edge), &
                    cell_vertices(:, next_edge), gauss_points(quadrature), &
                    point, tangent)
                contribution = point(1)*tangent(2) - &
                    point(2)*tangent(1)
                area = area + 0.5_dp*gauss_weights(quadrature)*contribution
            end do
        end do
    end function gauss_green_area

    pure subroutine bezier_point_and_tangent( &
            first, control, last, parameter, point, tangent)
        real(dp), intent(in) :: first(2), control(2), last(2), parameter
        real(dp), intent(out) :: point(2), tangent(2)
        real(dp) :: one_minus_parameter

        one_minus_parameter = 1.0_dp - parameter
        point = one_minus_parameter**2*first + &
            2.0_dp*one_minus_parameter*parameter*control + &
            parameter**2*last
        tangent = 2.0_dp*one_minus_parameter*(control - first) + &
            2.0_dp*parameter*(last - control)
    end subroutine bezier_point_and_tangent

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_fci_curved_polygon_area
