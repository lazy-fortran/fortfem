program fci_polygon_geometry
    ! FortSym-generated fixed-topology areas for arbitrary FCI polygons.
    ! The boundary-ordered pentagons are checked against independent
    ! shoelace and Gauss--Green oracles before their geometry is plotted.
    use fortfem_api, only: &
        compute_fci_cubic_curved_polygon_cell_areas_2d, &
        compute_fci_curved_polygon_cell_areas_2d, &
        compute_fci_polygon_cell_areas_2d
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: vertex_count = 5
    integer, parameter :: cell_count = 2
    integer, parameter :: edge_curve_points = 12
    real(dp), parameter :: gauss_points(3) = [ &
        0.1127016653792583_dp, 0.5_dp, 0.8872983346207417_dp]
    real(dp), parameter :: gauss_weights(3) = [ &
        5.0_dp/18.0_dp, 8.0_dp/18.0_dp, 5.0_dp/18.0_dp]
    character(*), parameter :: output_directory = &
        "output/example/fci_polygon_geometry"
    real(dp) :: vertices(2, vertex_count, cell_count)
    real(dp) :: edge_controls(2, vertex_count, cell_count)
    real(dp) :: cubic_controls(2, 2, vertex_count, cell_count)
    real(dp) :: areas(cell_count), expected(cell_count)
    real(dp) :: curved_areas(cell_count), curved_expected(cell_count)
    real(dp) :: cubic_areas(cell_count), cubic_expected(cell_count)
    real(dp) :: x_plot(vertex_count + 1), y_plot(vertex_count + 1)
    real(dp) :: x_curve(vertex_count*edge_curve_points + 1)
    real(dp) :: y_curve(vertex_count*edge_curve_points + 1)
    real(dp) :: x_cubic(vertex_count*edge_curve_points + 1)
    real(dp) :: y_cubic(vertex_count*edge_curve_points + 1)
    real(dp) :: bezier_point(2), bezier_tangent(2), parameter
    real(dp) :: cell_index(cell_count)
    integer :: cell, vertex, edge, sample, point_index, next_vertex, unit
    type(fortsparse_status_t) :: status

    call execute_command_line("mkdir -p "//output_directory)
    vertices(:, :, 1) = reshape([ &
        0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 1.0_dp, &
        1.6_dp, 2.0_dp, 0.0_dp, 1.4_dp], [2, vertex_count])
    vertices(:, :, 2) = reshape([ &
        4.0_dp, 0.0_dp, 6.0_dp, 0.3_dp, 6.4_dp, 1.7_dp, &
        5.0_dp, 2.3_dp, 3.9_dp, 1.1_dp], [2, vertex_count])
    edge_controls(:, :, 1) = reshape([ &
        1.0_dp, -0.25_dp, 2.7_dp, 0.5_dp, 2.4_dp, 1.7_dp, &
        0.7_dp, 2.1_dp, -0.2_dp, 0.6_dp], [2, vertex_count])
    edge_controls(:, :, 2) = reshape([ &
        5.0_dp, 0.0_dp, 6.6_dp, 1.0_dp, 5.6_dp, 2.4_dp, &
        4.2_dp, 1.8_dp, 3.8_dp, 0.4_dp], [2, vertex_count])
    cubic_controls(:, 1, :, :) = 0.0_dp
    cubic_controls(:, 2, :, :) = 0.0_dp
    cubic_controls(:, 1, :, 1) = edge_controls(:, :, 1)
    cubic_controls(:, 1, :, 2) = edge_controls(:, :, 2)
    cubic_controls(:, 2, :, 1) = reshape([ &
        1.5_dp, -0.35_dp, 3.0_dp, 0.8_dp, 2.7_dp, 1.9_dp, &
        0.4_dp, 2.2_dp, -0.1_dp, 1.0_dp], [2, vertex_count])
    cubic_controls(:, 2, :, 2) = reshape([ &
        5.5_dp, 0.1_dp, 6.7_dp, 1.3_dp, 5.3_dp, 2.6_dp, &
        4.0_dp, 2.0_dp, 3.7_dp, 0.7_dp], [2, vertex_count])

    call compute_fci_polygon_cell_areas_2d(vertices, areas, status)
    if (status%code /= 0) error stop "FCI polygon area evaluation failed"
    do cell = 1, cell_count
        expected(cell) = shoelace_area(vertices(:, :, cell))
        if (abs(areas(cell) - expected(cell)) > 3.0e-14_dp) then
            error stop "FCI polygon area oracle failed"
        end if
        cell_index(cell) = real(cell, dp)
    end do
    call compute_fci_curved_polygon_cell_areas_2d( &
        vertices, edge_controls, curved_areas, status)
    if (status%code /= 0) then
        error stop "FCI curved polygon area evaluation failed"
    end if
    do cell = 1, cell_count
        curved_expected(cell) = gauss_green_area( &
            vertices(:, :, cell), edge_controls(:, :, cell))
        if (abs(curved_areas(cell) - curved_expected(cell)) > 3.0e-12_dp) then
            error stop "FCI curved polygon area oracle failed"
        end if
    end do
    call compute_fci_cubic_curved_polygon_cell_areas_2d( &
        vertices, cubic_controls, cubic_areas, status)
    if (status%code /= 0) then
        error stop "FCI cubic curved polygon area evaluation failed"
    end if
    do cell = 1, cell_count
        cubic_expected(cell) = cubic_gauss_green_area( &
            vertices(:, :, cell), cubic_controls(:, :, :, cell))
        if (abs(cubic_areas(cell) - cubic_expected(cell)) > 3.0e-12_dp) then
            error stop "FCI cubic curved polygon area oracle failed"
        end if
    end do

    call figure(figsize=[9.0_dp, 6.0_dp])
    do cell = 1, cell_count
        do vertex = 1, vertex_count
            x_plot(vertex) = vertices(1, vertex, cell)
            y_plot(vertex) = vertices(2, vertex, cell)
        end do
        x_plot(vertex_count + 1) = x_plot(1)
        y_plot(vertex_count + 1) = y_plot(1)
        call plot(x_plot, y_plot, label="vertices "//integer_label(cell), &
            linestyle=":", marker="o")
        do edge = 1, vertex_count
            next_vertex = mod(edge, vertex_count) + 1
            do sample = 0, edge_curve_points
                point_index = (edge - 1)*edge_curve_points + sample + 1
                parameter = real(sample, dp)/real(edge_curve_points, dp)
                call cubic_bezier_point_and_tangent( &
                    vertices(:, edge, cell), cubic_controls(:, :, edge, cell), &
                    vertices(:, next_vertex, cell), parameter, bezier_point, &
                    bezier_tangent)
                x_cubic(point_index) = bezier_point(1)
                y_cubic(point_index) = bezier_point(2)
            end do
        end do
        call plot(x_cubic, y_cubic, label="cubic boundary "// &
            integer_label(cell), linestyle="-")
    end do
    call xlabel("poloidal x")
    call ylabel("poloidal y")
    call title("Cubic Bezier FCI polygon cells")
    call legend()
    call savefig(output_directory//"/fci_polygon_cells_2d.png")

    call figure(figsize=[9.0_dp, 6.0_dp])
    do cell = 1, cell_count
        do edge = 1, vertex_count
            next_vertex = mod(edge, vertex_count) + 1
            do sample = 0, edge_curve_points
                point_index = (edge - 1)*edge_curve_points + sample + 1
                parameter = real(sample, dp)/real(edge_curve_points, dp)
                call bezier_point_and_tangent( &
                    vertices(:, edge, cell), edge_controls(:, edge, cell), &
                    vertices(:, next_vertex, cell), parameter, bezier_point, &
                    bezier_tangent)
                x_curve(point_index) = bezier_point(1)
                y_curve(point_index) = bezier_point(2)
            end do
        end do
        call plot(x_curve, y_curve, label="curved pentagon "// &
            integer_label(cell))
        do edge = 1, vertex_count
            next_vertex = mod(edge, vertex_count) + 1
            do sample = 0, edge_curve_points
                point_index = (edge - 1)*edge_curve_points + sample + 1
                parameter = real(sample, dp)/real(edge_curve_points, dp)
                call cubic_bezier_point_and_tangent( &
                    vertices(:, edge, cell), cubic_controls(:, :, edge, cell), &
                    vertices(:, next_vertex, cell), parameter, bezier_point, &
                    bezier_tangent)
                x_cubic(point_index) = bezier_point(1)
                y_cubic(point_index) = bezier_point(2)
            end do
        end do
        call plot(x_cubic, y_cubic, label="cubic pentagon "// &
            integer_label(cell), linestyle="--")
    end do
    call xlabel("poloidal x")
    call ylabel("poloidal y")
    call title("Quadratic and cubic Bezier-edge FCI polygon cells")
    call legend()
    call savefig(output_directory//"/fci_curved_polygon_cells_2d.png")

    call figure(figsize=[8.0_dp, 5.0_dp])
    call plot(cell_index, areas, marker="o", label="generated area")
    call plot(cell_index, expected, marker="x", label="shoelace oracle")
    call plot(cell_index, curved_areas, marker="s", &
        label="quadratic Bezier-edge area")
    call plot(cell_index, cubic_areas, marker="^", linestyle="--", &
        label="cubic Bezier-edge area")
    call xlabel("cell index")
    call ylabel("cell area")
    call title("Generated FCI polygon plane-cell measures")
    call legend()
    call savefig(output_directory//"/fci_polygon_areas_1d.png")

    open (newunit=unit, file=output_directory//"/fci_polygon_areas.csv", &
        status="replace", action="write")
    write (unit, "(a)") "cell,generated_area,shoelace_area,quadratic_area,"// &
        "quadratic_gauss_green_area,cubic_area,cubic_gauss_green_area"
    do cell = 1, cell_count
        write (unit, "(i0,6(',',es24.16))") cell, areas(cell), expected(cell), &
            curved_areas(cell), curved_expected(cell), cubic_areas(cell), &
            cubic_expected(cell)
    end do
    close (unit)

contains

    pure real(dp) function shoelace_area(cell_vertices) result(area)
        real(dp), intent(in) :: cell_vertices(2, vertex_count)
        integer :: vertex, next_vertex

        area = 0.0_dp
        do vertex = 1, vertex_count
            next_vertex = mod(vertex, vertex_count) + 1
            area = area + cell_vertices(1, vertex)* &
                cell_vertices(2, next_vertex) - &
                cell_vertices(2, vertex)*cell_vertices(1, next_vertex)
        end do
        area = 0.5_dp*area
    end function shoelace_area

    pure real(dp) function gauss_green_area(cell_vertices, controls) result(area)
        real(dp), intent(in) :: cell_vertices(2, vertex_count)
        real(dp), intent(in) :: controls(2, vertex_count)
        real(dp) :: point(2), tangent(2), contribution
        integer :: edge, quadrature, next_edge

        area = 0.0_dp
        do edge = 1, vertex_count
            next_edge = mod(edge, vertex_count) + 1
            do quadrature = 1, 3
                call bezier_point_and_tangent( &
                    cell_vertices(:, edge), controls(:, edge), &
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

    pure real(dp) function cubic_gauss_green_area(cell_vertices, controls) &
            result(area)
        real(dp), intent(in) :: cell_vertices(2, vertex_count)
        real(dp), intent(in) :: controls(2, 2, vertex_count)
        real(dp) :: point(2), tangent(2), contribution
        integer :: edge, quadrature, next_edge

        area = 0.0_dp
        do edge = 1, vertex_count
            next_edge = mod(edge, vertex_count) + 1
            do quadrature = 1, 3
                call cubic_bezier_point_and_tangent( &
                    cell_vertices(:, edge), controls(:, :, edge), &
                    cell_vertices(:, next_edge), gauss_points(quadrature), &
                    point, tangent)
                contribution = point(1)*tangent(2) - &
                    point(2)*tangent(1)
                area = area + 0.5_dp*gauss_weights(quadrature)*contribution
            end do
        end do
    end function cubic_gauss_green_area

    pure subroutine cubic_bezier_point_and_tangent( &
            first, controls, last, parameter, point, tangent)
        real(dp), intent(in) :: first(2), controls(2, 2), last(2), parameter
        real(dp), intent(out) :: point(2), tangent(2)
        real(dp) :: one_minus_parameter

        one_minus_parameter = 1.0_dp - parameter
        point = one_minus_parameter**3*first + &
            3.0_dp*one_minus_parameter**2*parameter*controls(:, 1) + &
            3.0_dp*one_minus_parameter*parameter**2*controls(:, 2) + &
            parameter**3*last
        tangent = 3.0_dp*one_minus_parameter**2*(controls(:, 1) - first) + &
            6.0_dp*one_minus_parameter*parameter* &
            (controls(:, 2) - controls(:, 1)) + &
            3.0_dp*parameter**2*(last - controls(:, 2))
    end subroutine cubic_bezier_point_and_tangent

    function integer_label(value) result(label)
        integer, intent(in) :: value
        character(1) :: label

        write (label, "(i1)") value
    end function integer_label

end program fci_polygon_geometry
