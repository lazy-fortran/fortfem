program fci_quadrilateral_geometry
    ! FortSym-generated areas for unstructured FCI plane cells.
    use fortfem_api, only: compute_fci_quadrilateral_cell_areas_2d
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: cell_count = 3
    character(*), parameter :: output_directory = &
        "output/example/fci_quadrilateral_geometry"
    real(dp) :: vertices(2, 4, cell_count), areas(cell_count)
    real(dp) :: x_plot(5), y_plot(5), cell_index(cell_count)
    real(dp) :: expected_area
    integer :: cell, vertex, unit
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

    call compute_fci_quadrilateral_cell_areas_2d(vertices, areas, status)
    if (status%code /= 0) error stop "FCI quadrilateral area evaluation failed"
    do cell = 1, cell_count
        expected_area = shoelace_area(vertices(:, :, cell))
        if (abs(areas(cell) - expected_area) > 3.0e-14_dp) then
            error stop "FCI quadrilateral area oracle failed"
        end if
        cell_index(cell) = real(cell, dp)
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

    call figure(figsize=[8.0_dp, 5.0_dp])
    call plot(cell_index, areas, marker="o", label="signed shoelace area")
    call xlabel("cell index")
    call ylabel("cell area")
    call title("Generated FCI plane-cell measures")
    call legend()
    call savefig(output_directory//"/fci_quadrilateral_areas_1d.png")

    open (newunit=unit, file=output_directory//"/fci_quadrilateral_areas.csv", &
        status="replace", action="write")
    write (unit, "(a)") "cell,area"
    do cell = 1, cell_count
        write (unit, "(i0,',',es24.16)") cell, areas(cell)
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

    function integer_label(value) result(label)
        integer, intent(in) :: value
        character(1) :: label

        write (label, "(i1)") value
    end function integer_label

end program fci_quadrilateral_geometry
