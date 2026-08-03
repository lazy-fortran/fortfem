program test_fci_quadrilateral_area
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        compute_fci_quadrilateral_cell_areas_2d, &
        compute_fci_quadrilateral_cell_areas_2d_jvp, &
        compute_fci_quadrilateral_cell_areas_2d_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: cell_count = 2
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(2, 4, cell_count), vertices_dot(2, 4, cell_count)
    real(dp) :: vertices_bar(2, 4, cell_count)
    real(dp) :: vertices_plus(2, 4, cell_count), vertices_minus(2, 4, cell_count)
    real(dp) :: areas(cell_count), areas_dot(cell_count), areas_bar(cell_count)
    real(dp) :: areas_plus(cell_count), areas_minus(cell_count)
    real(dp) :: expected(cell_count), lhs, rhs
    real(dp) :: invalid_vertices(2, 4, cell_count)
    type(fortsparse_status_t) :: status
    integer :: cell, vertex
    logical :: all_passed

    all_passed = .true.
    vertices(:, :, 1) = reshape([ &
        0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 2.5_dp, 1.0_dp, 0.0_dp, 1.5_dp], [2, 4])
    vertices(:, :, 2) = reshape([ &
        1.0_dp, 0.2_dp, 2.4_dp, 0.4_dp, 2.0_dp, 1.6_dp, 0.8_dp, 1.3_dp], [2, 4])
    do cell = 1, cell_count
        do vertex = 1, 4
            vertices_dot(:, vertex, cell) = [ &
                0.01_dp*sin(real(3*vertex + cell, dp)), &
                -0.02_dp*cos(real(2*vertex + cell, dp))]
        end do
    end do
    vertices_bar = reshape([(0.03_dp*sin(real(5*vertex, dp)), &
        vertex=1, size(vertices_bar))], shape(vertices_bar))
    areas_bar = [0.7_dp, -1.1_dp]

    call compute_fci_quadrilateral_cell_areas_2d(vertices, areas, status)
    call record_condition(status%code == 0, &
        "FCI quadrilateral area accepts positively oriented cells")
    do cell = 1, cell_count
        expected(cell) = shoelace_area(vertices(:, :, cell))
    end do
    call record_condition(maxval(abs(areas - expected)) < 3.0e-14_dp, &
        "FCI quadrilateral area matches the independent shoelace oracle")
    call record_condition(all(areas > 0.0_dp), &
        "FCI quadrilateral areas are positive")

    call compute_fci_quadrilateral_cell_areas_2d_jvp( &
        vertices, vertices_dot, areas_dot, status)
    vertices_plus = vertices + step*vertices_dot
    vertices_minus = vertices - step*vertices_dot
    call compute_fci_quadrilateral_cell_areas_2d( &
        vertices_plus, areas_plus, status)
    call compute_fci_quadrilateral_cell_areas_2d( &
        vertices_minus, areas_minus, status)
    call record_condition(maxval(abs(areas_dot - &
        (areas_plus - areas_minus)/(2.0_dp*step))) < 3.0e-8_dp, &
        "FCI quadrilateral area JVP matches central differences")

    call compute_fci_quadrilateral_cell_areas_2d_vjp( &
        vertices, areas_bar, vertices_bar, status)
    lhs = dot_product(areas_bar, areas_dot)
    rhs = sum(vertices_bar*vertices_dot)
    call record_condition(abs(lhs - rhs) < &
        3.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "FCI quadrilateral area VJP satisfies the real dot-product identity")

    invalid_vertices = vertices
    invalid_vertices(:, :, 1) = vertices(:, [1, 4, 3, 2], 1)
    call compute_fci_quadrilateral_cell_areas_2d( &
        invalid_vertices, areas, status)
    call record_condition(status%code /= 0, &
        "FCI quadrilateral area rejects clockwise cells")
    invalid_vertices = vertices
    invalid_vertices(:, 4, 1) = invalid_vertices(:, 3, 1)
    call compute_fci_quadrilateral_cell_areas_2d( &
        invalid_vertices, areas, status)
    call record_condition(status%code /= 0, &
        "FCI quadrilateral area rejects degenerate cells")
    invalid_vertices = vertices
    invalid_vertices(:, :, 1) = reshape([ &
        0.0_dp, 0.0_dp, 3.0_dp, 3.0_dp, 0.0_dp, 2.0_dp, 2.0_dp, 0.0_dp], &
        [2, 4])
    call compute_fci_quadrilateral_cell_areas_2d( &
        invalid_vertices, areas, status)
    call record_condition(status%code /= 0, &
        "FCI quadrilateral area rejects self-intersecting cells")

    call check_summary("FCI quadrilateral area")
    if (.not. all_passed) error stop 1

contains

    pure real(dp) function shoelace_area(cell_vertices) result(area)
        real(dp), intent(in) :: cell_vertices(2, 4)
        integer :: vertex, next_vertex

        area = 0.0_dp
        do vertex = 1, 4
            next_vertex = mod(vertex, 4) + 1
            area = area + cell_vertices(1, vertex)*cell_vertices(2, next_vertex) - &
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

end program test_fci_quadrilateral_area
