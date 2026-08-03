program test_fci_septic_curved_polygon_area
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        compute_fci_septic_curved_polygon_cell_areas_2d, &
        compute_fci_septic_curved_polygon_cell_areas_2d_jvp, &
        compute_fci_septic_curved_polygon_cell_areas_2d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: vertex_count = 4
    integer, parameter :: cell_count = 2
    real(dp), parameter :: step = 1.0e-4_dp
    real(dp), parameter :: gauss_points(8) = [ &
        0.0198550717512319_dp, 0.1016667612931866_dp, &
        0.2372337950418355_dp, 0.4082826787521751_dp, &
        0.5917173212478249_dp, 0.7627662049581645_dp, &
        0.8983332387068134_dp, 0.9801449282487681_dp]
    real(dp), parameter :: gauss_weights(8) = [ &
        0.0506142681451881_dp, 0.1111905172266872_dp, &
        0.1568533229389436_dp, 0.1813418916891810_dp, &
        0.1813418916891810_dp, 0.1568533229389436_dp, &
        0.1111905172266872_dp, 0.0506142681451881_dp]
    real(dp) :: vertices(2, vertex_count, cell_count)
    real(dp) :: controls(2, 6, vertex_count, cell_count)
    real(dp) :: straight_controls(2, 6, vertex_count, cell_count)
    real(dp) :: bad_vertices(2, vertex_count, cell_count)
    real(dp) :: bad_controls(2, 6, vertex_count, cell_count)
    real(dp) :: vertices_dot(2, vertex_count, cell_count)
    real(dp) :: controls_dot(2, 6, vertex_count, cell_count)
    real(dp) :: vertices_bar(2, vertex_count, cell_count)
    real(dp) :: controls_bar(2, 6, vertex_count, cell_count)
    real(dp) :: vertices_plus(2, vertex_count, cell_count)
    real(dp) :: vertices_minus(2, vertex_count, cell_count)
    real(dp) :: controls_plus(2, 6, vertex_count, cell_count)
    real(dp) :: controls_minus(2, 6, vertex_count, cell_count)
    real(dp) :: areas(cell_count), areas_dot(cell_count), areas_bar(cell_count)
    real(dp) :: areas_plus(cell_count), areas_minus(cell_count)
    real(dp) :: expected(cell_count), straight(cell_count), chord(2), normal(2)
    real(dp) :: lhs, rhs, blend, curve_offset
    type(fortsparse_status_t) :: status
    integer :: cell, edge, control, next_edge
    logical :: all_passed

    all_passed = .true.
    vertices(:, :, 1) = reshape([ &
        0.0_dp, 0.0_dp, 2.0_dp, 0.1_dp, 2.5_dp, 1.6_dp, &
        0.1_dp, 2.0_dp], [2, vertex_count])
    vertices(:, :, 2) = reshape([ &
        3.7_dp, -0.2_dp, 6.2_dp, 0.2_dp, 5.9_dp, 2.3_dp, &
        3.5_dp, 1.4_dp], [2, vertex_count])
    do cell = 1, cell_count
        do edge = 1, vertex_count
            next_edge = mod(edge, vertex_count) + 1
            chord = vertices(:, next_edge, cell) - vertices(:, edge, cell)
            normal = [-chord(2), chord(1)]
            vertices_dot(:, edge, cell) = [ &
                0.013_dp*sin(real(2*edge + cell, dp)), &
                -0.017_dp*cos(real(3*edge + cell, dp))]
            do control = 1, 6
                blend = real(control, dp)/7.0_dp
                straight_controls(:, control, edge, cell) = &
                    vertices(:, edge, cell) + blend*chord
                curve_offset = 0.025_dp*real((-1)**(edge + cell), dp)* &
                    sin(acos(-1.0_dp)*blend)
                controls(:, control, edge, cell) = &
                    straight_controls(:, control, edge, cell) + &
                    curve_offset*normal
                controls_dot(:, control, edge, cell) = [ &
                    -0.011_dp*cos(real(edge + control + cell, dp)), &
                    0.014_dp*sin(real(2*edge + control + cell, dp))]
            end do
        end do
    end do
    areas_bar = [0.8_dp, -1.3_dp]

    call compute_fci_septic_curved_polygon_cell_areas_2d( &
        vertices, controls, areas, status)
    call record_condition(status%code == 0, &
        "FCI septic curved polygon area accepts curved cells")
    do cell = 1, cell_count
        expected(cell) = gauss_green_area( &
            vertices(:, :, cell), controls(:, :, :, cell))
        straight(cell) = shoelace_area(vertices(:, :, cell))
    end do
    call record_condition(maxval(abs(areas - expected)) < 8.0e-12_dp, &
        "FCI septic curved polygon area matches Gauss-Green oracle")
    call record_condition(maxval(abs(areas - straight)) > 1.0e-3_dp, &
        "FCI septic curved polygon area resolves edge curvature")
    call record_condition(all(areas > 0.0_dp), &
        "FCI septic curved polygon areas are positive")

    call compute_fci_septic_curved_polygon_cell_areas_2d_jvp( &
        vertices, controls, vertices_dot, controls_dot, areas_dot, status)
    vertices_plus = vertices + step*vertices_dot
    vertices_minus = vertices - step*vertices_dot
    controls_plus = controls + step*controls_dot
    controls_minus = controls - step*controls_dot
    call compute_fci_septic_curved_polygon_cell_areas_2d( &
        vertices_plus, controls_plus, areas_plus, status)
    call compute_fci_septic_curved_polygon_cell_areas_2d( &
        vertices_minus, controls_minus, areas_minus, status)
    call record_condition(maxval(abs(areas_dot - &
        (areas_plus - areas_minus)/(2.0_dp*step))) < 1.0e-7_dp, &
        "FCI septic curved polygon area JVP matches central differences")

    call compute_fci_septic_curved_polygon_cell_areas_2d_vjp( &
        vertices, controls, areas_bar, vertices_bar, controls_bar, status)
    lhs = dot_product(areas_bar, areas_dot)
    rhs = sum(vertices_bar*vertices_dot) + sum(controls_bar*controls_dot)
    call record_condition(abs(lhs - rhs) < &
        1.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "FCI septic curved polygon area VJP satisfies the real dot-product identity")

    call compute_fci_septic_curved_polygon_cell_areas_2d( &
        vertices, straight_controls, areas_plus, status)
    call record_condition(maxval(abs(areas_plus - straight)) < 8.0e-12_dp, &
        "FCI septic curved area recovers straight polygon area")

    bad_vertices = vertices
    bad_controls = straight_controls
    do edge = 1, vertex_count
        bad_vertices(:, edge, 1) = [real(edge - 1, dp), 0.0_dp]
    end do
    do edge = 1, vertex_count
        next_edge = mod(edge, vertex_count) + 1
        do control = 1, 6
            blend = real(control, dp)/7.0_dp
            bad_controls(:, control, edge, 1) = bad_vertices(:, edge, 1) + &
                blend*(bad_vertices(:, next_edge, 1) - bad_vertices(:, edge, 1))
        end do
    end do
    call compute_fci_septic_curved_polygon_cell_areas_2d( &
        bad_vertices, bad_controls, areas_plus, status)
    call record_condition(status%code /= 0 .and. all(areas_plus == 0.0_dp), &
        "FCI septic curved polygon area rejects degenerate cells")

    call check_summary("FCI septic curved polygon area")
    if (.not. all_passed) error stop 1

contains

    pure real(dp) function shoelace_area(cell_vertices) result(area)
        real(dp), intent(in) :: cell_vertices(2, vertex_count)
        integer :: edge, next_edge

        area = 0.0_dp
        do edge = 1, vertex_count
            next_edge = mod(edge, vertex_count) + 1
            area = area + cell_vertices(1, edge)*cell_vertices(2, next_edge) - &
                cell_vertices(2, edge)*cell_vertices(1, next_edge)
        end do
        area = 0.5_dp*area
    end function shoelace_area

    pure real(dp) function gauss_green_area(cell_vertices, edge_controls) &
            result(area)
        real(dp), intent(in) :: cell_vertices(2, vertex_count)
        real(dp), intent(in) :: edge_controls(2, 6, vertex_count)
        real(dp) :: point(2), tangent(2)
        integer :: edge, quadrature, next_edge

        area = 0.0_dp
        do edge = 1, vertex_count
            next_edge = mod(edge, vertex_count) + 1
            do quadrature = 1, size(gauss_points)
                call septic_bezier_point_and_tangent( &
                    cell_vertices(:, edge), edge_controls(:, :, edge), &
                    cell_vertices(:, next_edge), gauss_points(quadrature), &
                    point, tangent)
                area = area + 0.5_dp*gauss_weights(quadrature)* &
                    (point(1)*tangent(2) - point(2)*tangent(1))
            end do
        end do
    end function gauss_green_area

    pure subroutine septic_bezier_point_and_tangent( &
            first, controls, last, parameter, point, tangent)
        real(dp), intent(in) :: first(2), controls(2, 6), last(2), parameter
        real(dp), intent(out) :: point(2), tangent(2)
        real(dp) :: bernstein(0:7), derivative_bernstein(0:6)
        real(dp) :: control_points(2, 0:7), one_minus_parameter
        integer :: index

        one_minus_parameter = 1.0_dp - parameter
        control_points(:, 0) = first
        control_points(:, 1:6) = controls
        control_points(:, 7) = last
        bernstein = [ &
            one_minus_parameter**7, &
            7.0_dp*one_minus_parameter**6*parameter, &
            21.0_dp*one_minus_parameter**5*parameter**2, &
            35.0_dp*one_minus_parameter**4*parameter**3, &
            35.0_dp*one_minus_parameter**3*parameter**4, &
            21.0_dp*one_minus_parameter**2*parameter**5, &
            7.0_dp*one_minus_parameter*parameter**6, parameter**7]
        derivative_bernstein = [ &
            one_minus_parameter**6, &
            6.0_dp*one_minus_parameter**5*parameter, &
            15.0_dp*one_minus_parameter**4*parameter**2, &
            20.0_dp*one_minus_parameter**3*parameter**3, &
            15.0_dp*one_minus_parameter**2*parameter**4, &
            6.0_dp*one_minus_parameter*parameter**5, parameter**6]
        point = 0.0_dp
        tangent = 0.0_dp
        do index = 0, 7
            point = point + bernstein(index)*control_points(:, index)
        end do
        do index = 0, 6
            tangent = tangent + 7.0_dp*derivative_bernstein(index)* &
                (control_points(:, index + 1) - control_points(:, index))
        end do
    end subroutine septic_bezier_point_and_tangent

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_fci_septic_curved_polygon_area
