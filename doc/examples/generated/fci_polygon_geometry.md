---
title: fci_polygon_geometry Example
---

# fci_polygon_geometry Example

---
title: FCI polygon geometry
---

This fixture computes fixed-topology areas for boundary-ordered pentagonal FCI
plane cells with a FortSym-generated edge contribution.  It compares the
generated values with an independent shoelace oracle and writes a 2D polygon
boundary plot, a 1D area comparison, and CSV data for the gallery.

## Usage

```bash
fpm run --example fci_polygon_geometry
```

## Source Code

```fortran
program fci_polygon_geometry
    ! FortSym-generated fixed-topology areas for arbitrary FCI polygons.
    ! The boundary-ordered pentagons are checked against an independent
    ! shoelace oracle before their geometry is plotted for the gallery.
    use fortfem_api, only: compute_fci_polygon_cell_areas_2d
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: vertex_count = 5
    integer, parameter :: cell_count = 2
    character(*), parameter :: output_directory = &
        "output/example/fci_polygon_geometry"
    real(dp) :: vertices(2, vertex_count, cell_count)
    real(dp) :: areas(cell_count), expected(cell_count)
    real(dp) :: x_plot(vertex_count + 1), y_plot(vertex_count + 1)
    real(dp) :: cell_index(cell_count)
    integer :: cell, vertex, unit
    type(fortsparse_status_t) :: status

    call execute_command_line("mkdir -p "//output_directory)
    vertices(:, :, 1) = reshape([ &
        0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 1.0_dp, &
        1.6_dp, 2.0_dp, 0.0_dp, 1.4_dp], [2, vertex_count])
    vertices(:, :, 2) = reshape([ &
        4.0_dp, 0.0_dp, 6.0_dp, 0.3_dp, 6.4_dp, 1.7_dp, &
        5.0_dp, 2.3_dp, 3.9_dp, 1.1_dp], [2, vertex_count])

    call compute_fci_polygon_cell_areas_2d(vertices, areas, status)
    if (status%code /= 0) error stop "FCI polygon area evaluation failed"
    do cell = 1, cell_count
        expected(cell) = shoelace_area(vertices(:, :, cell))
        if (abs(areas(cell) - expected(cell)) > 3.0e-14_dp) then
            error stop "FCI polygon area oracle failed"
        end if
        cell_index(cell) = real(cell, dp)
    end do

    call figure(figsize=[9.0_dp, 6.0_dp])
    do cell = 1, cell_count
        do vertex = 1, vertex_count
            x_plot(vertex) = vertices(1, vertex, cell)
            y_plot(vertex) = vertices(2, vertex, cell)
        end do
        x_plot(vertex_count + 1) = x_plot(1)
        y_plot(vertex_count + 1) = y_plot(1)
        call plot(x_plot, y_plot, label="pentagon "//integer_label(cell))
    end do
    call xlabel("poloidal x")
    call ylabel("poloidal y")
    call title("Boundary-ordered FCI polygon cells")
    call legend()
    call savefig(output_directory//"/fci_polygon_cells_2d.png")

    call figure(figsize=[8.0_dp, 5.0_dp])
    call plot(cell_index, areas, marker="o", label="generated area")
    call plot(cell_index, expected, marker="x", label="shoelace oracle")
    call xlabel("cell index")
    call ylabel("cell area")
    call title("Generated FCI polygon plane-cell measures")
    call legend()
    call savefig(output_directory//"/fci_polygon_areas_1d.png")

    open (newunit=unit, file=output_directory//"/fci_polygon_areas.csv", &
        status="replace", action="write")
    write (unit, "(a)") "cell,generated_area,shoelace_area"
    do cell = 1, cell_count
        write (unit, "(i0,2(',',es24.16))") cell, areas(cell), expected(cell)
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

    function integer_label(value) result(label)
        integer, intent(in) :: value
        character(1) :: label

        write (label, "(i1)") value
    end function integer_label

end program fci_polygon_geometry
```

## Generated Plots

### cover.svg

![cover.svg](../../media/examples/fci_polygon_geometry/cover.svg)

---

[← Back to all examples](../index.html)
