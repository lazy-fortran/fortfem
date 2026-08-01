---
title: fci_quadrilateral_geometry Example
---

# fci_quadrilateral_geometry Example

---
title: FCI quadrilateral geometry
---

This small geometry fixture computes positive, fixed-topology areas for three
unstructured counter-clockwise quadrilateral FCI plane cells, including
quadratic Bezier-edge curvature.  It compares the generated values with
independent shoelace and Gauss--Green oracles and writes straight and curved
2D cell plots, an area summary plot, and CSV data for the gallery.

## Usage

```bash
fpm run --example fci_quadrilateral_geometry
```

## Source Code

```fortran
program fci_quadrilateral_geometry
    ! FortSym-generated straight and curved areas for FCI plane cells.
    use fortfem_api, only: &
        compute_fci_curved_quadrilateral_cell_areas_2d, &
        compute_fci_quadrilateral_cell_areas_2d
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: cell_count = 3
    integer, parameter :: edge_curve_points = 12
    real(dp), parameter :: gauss_points(3) = [ &
        0.1127016653792583_dp, 0.5_dp, 0.8872983346207417_dp]
    real(dp), parameter :: gauss_weights(3) = [ &
        5.0_dp/18.0_dp, 8.0_dp/18.0_dp, 5.0_dp/18.0_dp]
    character(*), parameter :: output_directory = &
        "output/example/fci_quadrilateral_geometry"
    real(dp) :: vertices(2, 4, cell_count), areas(cell_count)
    real(dp) :: edge_controls(2, 4, cell_count)
    real(dp) :: curved_areas(cell_count)
    real(dp) :: x_plot(5), y_plot(5), cell_index(cell_count)
    real(dp) :: x_curve(4*edge_curve_points + 1)
    real(dp) :: y_curve(4*edge_curve_points + 1)
    real(dp) :: bezier_point(2), bezier_tangent(2)
    real(dp) :: expected_area, parameter
    integer :: cell, edge, vertex, sample, point_index, next_edge, unit
    type(fortsparse_status_t) :: status

    call execute_command_line("mkdir -p "//output_directory)
    vertices(:, :, 1) = reshape([ &
        0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 2.5_dp, 1.0_dp, 0.0_dp, 1.5_dp], &
        [2, 4])
    vertices(:, :, 2) = reshape([ &
        2.8_dp, 0.2_dp, 4.4_dp, 0.4_dp, 4.0_dp, 1.8_dp, 2.5_dp, 1.3_dp], &
        [2, 4])
    vertices(:, :, 3) = reshape([ &
        0.4_dp, 2.1_dp, 1.8_dp, 2.0_dp, 2.3_dp, 3.3_dp, 0.1_dp, 3.0_dp], &
        [2, 4])
    edge_controls(:, :, 1) = reshape([ &
        1.0_dp, -0.25_dp, 2.35_dp, 0.45_dp, 1.25_dp, 1.55_dp, &
        -0.2_dp, 0.75_dp], [2, 4])
    edge_controls(:, :, 2) = reshape([ &
        3.6_dp, 0.0_dp, 4.3_dp, 1.1_dp, 3.2_dp, 1.9_dp, &
        2.6_dp, 0.75_dp], [2, 4])
    edge_controls(:, :, 3) = reshape([ &
        1.1_dp, 1.9_dp, 2.25_dp, 2.65_dp, 1.1_dp, 3.35_dp, &
        0.15_dp, 2.55_dp], [2, 4])

    call compute_fci_quadrilateral_cell_areas_2d(vertices, areas, status)
    if (status%code /= 0) error stop "FCI quadrilateral area evaluation failed"
    do cell = 1, cell_count
        expected_area = shoelace_area(vertices(:, :, cell))
        if (abs(areas(cell) - expected_area) > 3.0e-14_dp) then
            error stop "FCI quadrilateral area oracle failed"
        end if
        cell_index(cell) = real(cell, dp)
    end do
    call compute_fci_curved_quadrilateral_cell_areas_2d( &
        vertices, edge_controls, curved_areas, status)
    if (status%code /= 0) then
        error stop "FCI curved quadrilateral area evaluation failed"
    end if
    do cell = 1, cell_count
        expected_area = curved_green_area(vertices(:, :, cell), &
            edge_controls(:, :, cell))
        if (abs(curved_areas(cell) - expected_area) > 3.0e-12_dp) then
            error stop "FCI curved quadrilateral area oracle failed"
        end if
    end do

    call figure(figsize=[9.0_dp, 6.0_dp])
    do cell = 1, cell_count
        do vertex = 1, 4
            x_plot(vertex) = vertices(1, vertex, cell)
            y_plot(vertex) = vertices(2, vertex, cell)
        end do
        x_plot(5) = x_plot(1)
        y_plot(5) = y_plot(1)
        call plot(x_plot, y_plot, label="cell "//integer_label(cell))
    end do
    call xlabel("poloidal x")
    call ylabel("poloidal y")
    call title("Unstructured FCI quadrilateral plane cells")
    call legend()
    call savefig(output_directory//"/fci_quadrilateral_cells_2d.png")

    call figure(figsize=[9.0_dp, 6.0_dp])
    do cell = 1, cell_count
        do edge = 1, 4
            next_edge = mod(edge, 4) + 1
            do sample = 0, edge_curve_points
                point_index = (edge - 1)*edge_curve_points + sample + 1
                parameter = real(sample, dp)/real(edge_curve_points, dp)
                call bezier_point_and_tangent( &
                    vertices(:, edge, cell), edge_controls(:, edge, cell), &
                    vertices(:, next_edge, cell), parameter, bezier_point, &
                    bezier_tangent)
                x_curve(point_index) = bezier_point(1)
                y_curve(point_index) = bezier_point(2)
            end do
        end do
        call plot(x_curve, y_curve, label="curved cell "//integer_label(cell))
    end do
    call xlabel("poloidal x")
    call ylabel("poloidal y")
    call title("Quadratic Bezier-edge FCI plane cells")
    call legend()
    call savefig(output_directory//"/fci_curved_quadrilateral_cells_2d.png")

    call figure(figsize=[8.0_dp, 5.0_dp])
    call plot(cell_index, areas, marker="o", label="signed shoelace area")
    call plot(cell_index, curved_areas, marker="s", &
        label="quadratic Bezier-edge area")
    call xlabel("cell index")
    call ylabel("cell area")
    call title("Generated FCI plane-cell measures")
    call legend()
    call savefig(output_directory//"/fci_quadrilateral_areas_1d.png")

    open (newunit=unit, file=output_directory//"/fci_quadrilateral_areas.csv", &
        status="replace", action="write")
    write (unit, "(a)") "cell,straight_area,curved_area"
    do cell = 1, cell_count
        write (unit, "(i0,2(',',es24.16))") cell, areas(cell), curved_areas(cell)
    end do
    close (unit)

contains

    pure real(dp) function shoelace_area(cell_vertices) result(area)
        real(dp), intent(in) :: cell_vertices(2, 4)
        integer :: vertex, next_vertex

        area = 0.0_dp
        do vertex = 1, 4
            next_vertex = mod(vertex, 4) + 1
            area = area + cell_vertices(1, vertex)* &
                cell_vertices(2, next_vertex) - &
                cell_vertices(2, vertex)*cell_vertices(1, next_vertex)
        end do
        area = 0.5_dp*area
    end function shoelace_area

    pure real(dp) function curved_green_area(cell_vertices, controls) result(area)
        real(dp), intent(in) :: cell_vertices(2, 4), controls(2, 4)
        real(dp) :: point(2), tangent(2), contribution
        integer :: edge, quadrature, next_edge

        area = 0.0_dp
        do edge = 1, 4
            next_edge = mod(edge, 4) + 1
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
    end function curved_green_area

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

    function integer_label(value) result(label)
        integer, intent(in) :: value
        character(1) :: label

        write (label, "(i1)") value
    end function integer_label

end program fci_quadrilateral_geometry
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/fci_quadrilateral_geometry/primary.png)

### fci_curved_quadrilateral_cells_2d.png

![fci_curved_quadrilateral_cells_2d.png](../../media/examples/fci_quadrilateral_geometry/fci_curved_quadrilateral_cells_2d.png)

### fci_quadrilateral_areas_1d.png

![fci_quadrilateral_areas_1d.png](../../media/examples/fci_quadrilateral_geometry/fci_quadrilateral_areas_1d.png)

### fci_quadrilateral_cells_2d.png

![fci_quadrilateral_cells_2d.png](../../media/examples/fci_quadrilateral_geometry/fci_quadrilateral_cells_2d.png)

---

[← Back to all examples](../index.html)
