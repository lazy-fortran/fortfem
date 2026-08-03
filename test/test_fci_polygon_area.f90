program test_fci_polygon_area
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        compute_fci_polygon_cell_areas_2d, &
        compute_fci_polygon_cell_areas_2d_jvp, &
        compute_fci_polygon_cell_areas_2d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: vertex_count = 5
    integer, parameter :: cell_count = 2
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(2, vertex_count, cell_count)
    real(dp) :: vertices_dot(2, vertex_count, cell_count)
    real(dp) :: vertices_bar(2, vertex_count, cell_count)
    real(dp) :: vertices_plus(2, vertex_count, cell_count)
    real(dp) :: vertices_minus(2, vertex_count, cell_count)
    real(dp) :: areas(cell_count), areas_dot(cell_count), areas_bar(cell_count)
    real(dp) :: areas_plus(cell_count), areas_minus(cell_count)
    real(dp) :: expected(cell_count), lhs, rhs
    real(dp) :: invalid_vertices(2, vertex_count, cell_count)
    type(fortsparse_status_t) :: status
    integer :: cell, vertex
    logical :: all_passed

    all_passed = .true.
    vertices(:, :, 1) = reshape([ &
        0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 1.0_dp, &
        1.6_dp, 2.0_dp, 0.0_dp, 1.4_dp], [2, vertex_count])
    vertices(:, :, 2) = reshape([ &
        4.0_dp, 0.0_dp, 6.0_dp, 0.3_dp, 6.4_dp, 1.7_dp, &
        5.0_dp, 2.3_dp, 3.9_dp, 1.1_dp], [2, vertex_count])
    do cell = 1, cell_count
        do vertex = 1, vertex_count
            vertices_dot(:, vertex, cell) = [ &
                0.01_dp*sin(real(3*vertex + cell, dp)), &
                -0.02_dp*cos(real(2*vertex + cell, dp))]
        end do
    end do
    vertices_bar = reshape([(0.03_dp*sin(real(5*vertex, dp)), &
        vertex=1, size(vertices_bar))], shape(vertices_bar))
    areas_bar = [0.7_dp, -1.1_dp]

    call compute_fci_polygon_cell_areas_2d(vertices, areas, status)
    call record_condition(status%code == 0, &
        "FCI polygon area accepts positively oriented pentagons")
    do cell = 1, cell_count
        expected(cell) = shoelace_area(vertices(:, :, cell))
    end do
    call record_condition(maxval(abs(areas - expected)) < 3.0e-14_dp, &
        "FCI polygon area matches the independent shoelace oracle")
    call record_condition(all(areas > 0.0_dp), &
        "FCI polygon areas are positive")

    call compute_fci_polygon_cell_areas_2d_jvp( &
        vertices, vertices_dot, areas_dot, status)
    vertices_plus = vertices + step*vertices_dot
    vertices_minus = vertices - step*vertices_dot
    call compute_fci_polygon_cell_areas_2d( &
        vertices_plus, areas_plus, status)
    call compute_fci_polygon_cell_areas_2d( &
        vertices_minus, areas_minus, status)
    call record_condition(maxval(abs(areas_dot - &
        (areas_plus - areas_minus)/(2.0_dp*step))) < 3.0e-8_dp, &
        "FCI polygon area JVP matches central differences")

    call compute_fci_polygon_cell_areas_2d_vjp( &
        vertices, areas_bar, vertices_bar, status)
    lhs = dot_product(areas_bar, areas_dot)
    rhs = sum(vertices_bar*vertices_dot)
    call record_condition(abs(lhs - rhs) < &
        3.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "FCI polygon area VJP satisfies the real dot-product identity")

    invalid_vertices = vertices
    invalid_vertices(:, :, 1) = vertices(:, [1, 5, 4, 3, 2], 1)
    call compute_fci_polygon_cell_areas_2d( &
        invalid_vertices, areas, status)
    call record_condition(status%code /= 0, &
        "FCI polygon area rejects clockwise cells")

    invalid_vertices = vertices
    invalid_vertices(:, 5, 1) = invalid_vertices(:, 4, 1)
    call compute_fci_polygon_cell_areas_2d( &
        invalid_vertices, areas, status)
    call record_condition(status%code /= 0, &
        "FCI polygon area rejects repeated vertices")

    invalid_vertices = vertices
    invalid_vertices(:, :, 1) = reshape([ &
        0.0_dp, 0.0_dp, 3.0_dp, 3.0_dp, 0.0_dp, 3.0_dp, &
        3.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [2, vertex_count])
    call compute_fci_polygon_cell_areas_2d( &
        invalid_vertices, areas, status)
    call record_condition(status%code /= 0, &
        "FCI polygon area rejects self-intersecting cells")

    call check_summary("FCI polygon area")
    if (.not. all_passed) error stop 1

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

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_fci_polygon_area
