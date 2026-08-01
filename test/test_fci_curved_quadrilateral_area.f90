program test_fci_curved_quadrilateral_area
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        compute_fci_curved_quadrilateral_cell_areas_2d, &
        compute_fci_curved_quadrilateral_cell_areas_2d_jvp, &
        compute_fci_curved_quadrilateral_cell_areas_2d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: cell_count = 2
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp), parameter :: gauss_points(3) = [ &
        0.1127016653792583_dp, 0.5_dp, 0.8872983346207417_dp]
    real(dp), parameter :: gauss_weights(3) = [ &
        5.0_dp/18.0_dp, 8.0_dp/18.0_dp, 5.0_dp/18.0_dp]
    real(dp) :: vertices(2, 4, cell_count), controls(2, 4, cell_count)
    real(dp) :: vertices_dot(2, 4, cell_count), controls_dot(2, 4, cell_count)
    real(dp) :: vertices_bar(2, 4, cell_count), controls_bar(2, 4, cell_count)
    real(dp) :: vertices_plus(2, 4, cell_count), vertices_minus(2, 4, cell_count)
    real(dp) :: controls_plus(2, 4, cell_count), controls_minus(2, 4, cell_count)
    real(dp) :: areas(cell_count), areas_dot(cell_count), areas_bar(cell_count)
    real(dp) :: areas_plus(cell_count), areas_minus(cell_count)
    real(dp) :: expected(cell_count), straight(cell_count), lhs, rhs
    real(dp) :: midpoint_controls(2, 4, cell_count)
    type(fortsparse_status_t) :: status
    integer :: cell, edge, next_edge
    logical :: all_passed

    all_passed = .true.
    vertices(:, :, 1) = reshape([ &
        0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 2.5_dp, 1.0_dp, 0.0_dp, 1.5_dp], &
        [2, 4])
    vertices(:, :, 2) = reshape([ &
        1.0_dp, 2.0_dp, 2.6_dp, 2.1_dp, 2.2_dp, 3.6_dp, 0.6_dp, 3.3_dp], &
        [2, 4])
    controls(:, :, 1) = reshape([ &
        1.0_dp, -0.25_dp, 2.35_dp, 0.45_dp, 1.25_dp, 1.55_dp, &
        -0.2_dp, 0.75_dp], [2, 4])
    controls(:, :, 2) = reshape([ &
        1.8_dp, 1.75_dp, 2.55_dp, 2.85_dp, 1.35_dp, 3.55_dp, &
        0.35_dp, 2.55_dp], [2, 4])
    do cell = 1, cell_count
        do edge = 1, 4
            next_edge = mod(edge, 4) + 1
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

    call compute_fci_curved_quadrilateral_cell_areas_2d( &
        vertices, controls, areas, status)
    call record_condition(status%code == 0, &
        "FCI curved quadrilateral area accepts curved cells")
    do cell = 1, cell_count
        expected(cell) = gauss_green_area(vertices(:, :, cell), &
            controls(:, :, cell))
        straight(cell) = polygon_area(vertices(:, :, cell))
    end do
    call record_condition(maxval(abs(areas - expected)) < 2.0e-12_dp, &
        "FCI curved area matches independent Gauss Green oracle")
    call record_condition(maxval(abs(areas - straight)) > 1.0e-3_dp, &
        "FCI curved area reflects nonzero edge curvature")
    call record_condition(all(areas > 0.0_dp), &
        "FCI curved quadrilateral areas are positive")

    call compute_fci_curved_quadrilateral_cell_areas_2d_jvp( &
        vertices, controls, vertices_dot, controls_dot, areas_dot, status)
    vertices_plus = vertices + step*vertices_dot
    vertices_minus = vertices - step*vertices_dot
    controls_plus = controls + step*controls_dot
    controls_minus = controls - step*controls_dot
    call compute_fci_curved_quadrilateral_cell_areas_2d( &
        vertices_plus, controls_plus, areas_plus, status)
    call compute_fci_curved_quadrilateral_cell_areas_2d( &
        vertices_minus, controls_minus, areas_minus, status)
    call record_condition(maxval(abs(areas_dot - &
        (areas_plus - areas_minus)/(2.0_dp*step))) < 4.0e-8_dp, &
        "FCI curved area JVP matches central differences")

    call compute_fci_curved_quadrilateral_cell_areas_2d_vjp( &
        vertices, controls, areas_bar, vertices_bar, controls_bar, status)
    lhs = dot_product(areas_bar, areas_dot)
    rhs = sum(vertices_bar*vertices_dot) + sum(controls_bar*controls_dot)
    call record_condition(abs(lhs - rhs) < &
        5.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "FCI curved area VJP satisfies the real dot-product identity")

    call compute_fci_curved_quadrilateral_cell_areas_2d( &
        vertices, midpoint_controls, straight, status)
    call record_condition(maxval(abs(straight - &
        [(polygon_area(vertices(:, :, cell)), cell=1, cell_count)])) < &
        2.0e-12_dp, "FCI curved area recovers polygon area for straight edges")

    vertices(:, 4, 1) = vertices(:, 3, 1)
    call compute_fci_curved_quadrilateral_cell_areas_2d( &
        vertices, controls, areas, status)
    call record_condition(status%code /= 0, &
        "FCI curved area rejects invalid curved-cell topology")

    call check_summary("FCI curved quadrilateral area")
    if (.not. all_passed) error stop 1

contains

    pure real(dp) function polygon_area(cell_vertices) result(area)
        real(dp), intent(in) :: cell_vertices(2, 4)
        integer :: edge, next_edge

        area = 0.0_dp
        do edge = 1, 4
            next_edge = mod(edge, 4) + 1
            area = area + cell_vertices(1, edge)* &
                cell_vertices(2, next_edge) - &
                cell_vertices(2, edge)*cell_vertices(1, next_edge)
        end do
        area = 0.5_dp*area
    end function polygon_area

    pure real(dp) function gauss_green_area(cell_vertices, edge_controls) &
            result(area)
        real(dp), intent(in) :: cell_vertices(2, 4), edge_controls(2, 4)
        real(dp) :: point(2), tangent(2), contribution
        integer :: edge, quadrature

        area = 0.0_dp
        do edge = 1, 4
            do quadrature = 1, 3
                call bezier_point_and_tangent( &
                    cell_vertices(:, edge), edge_controls(:, edge), &
                    cell_vertices(:, mod(edge, 4) + 1), &
                    gauss_points(quadrature), point, tangent)
                contribution = point(1)*tangent(2) - point(2)*tangent(1)
                area = area + 0.5_dp*gauss_weights(quadrature)*contribution
            end do
        end do
    end function gauss_green_area

    pure subroutine bezier_point_and_tangent(first, control, last, parameter, &
            point, tangent)
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

end program test_fci_curved_quadrilateral_area
