program test_region_interface_graph
    use check, only: check_condition, check_summary
    use fortfem_api, only: initialize_region_interface_graph, &
        region_interface_graph_components, region_interface_graph_incidence, &
        region_interface_graph_t, validate_region_interface_graph
    implicit none

    type(region_interface_graph_t) :: chain, disconnected, periodic, invalid
    integer, allocatable :: plus_region(:), minus_region(:)
    integer, allocatable :: incidence(:, :), components(:)
    integer :: component_count, status

    allocate(plus_region(2), minus_region(2))
    plus_region = [2, 3]
    minus_region = [1, 2]
    call initialize_region_interface_graph( &
        chain, 3, plus_region, minus_region, status)
    call check_condition(status == 0, "oriented region chain initializes")
    call validate_region_interface_graph(chain, status)
    call check_condition(status == 0, "oriented region chain validates")
    call region_interface_graph_incidence(chain, incidence, status)
    call check_condition(status == 0 .and. all(incidence == reshape([ &
        -1, 1, 0, 0, -1, 1], [3, 2])), &
        "region incidence has the independent plus-minus oracle")
    call region_interface_graph_components( &
        chain, components, component_count, status)
    call check_condition(status == 0 .and. component_count == 1 .and. &
        all(components == [1, 1, 1]), &
        "connected region chain has one component")

    plus_region = [1, 2]
    minus_region = [2, 3]
    call initialize_region_interface_graph( &
        chain, 3, plus_region, minus_region, status)
    call region_interface_graph_incidence(chain, incidence, status)
    call check_condition(status == 0 .and. all(incidence == reshape([ &
        1, -1, 0, 0, 1, -1], [3, 2])), &
        "reversing each interface reverses its incidence column")

    deallocate(plus_region, minus_region)
    allocate(plus_region(2), minus_region(2))
    plus_region = [2, 4]
    minus_region = [1, 3]
    call initialize_region_interface_graph( &
        disconnected, 4, plus_region, minus_region, status)
    call region_interface_graph_components( &
        disconnected, components, component_count, status)
    call check_condition(status == 0 .and. component_count == 2 .and. &
        all(components == [1, 1, 2, 2]), &
        "disconnected region graph preserves two components")

    deallocate(plus_region, minus_region)
    allocate(plus_region(1), minus_region(1))
    plus_region = 1
    minus_region = 1
    call initialize_region_interface_graph( &
        periodic, 1, plus_region, minus_region, status)
    call region_interface_graph_incidence(periodic, incidence, status)
    call check_condition(status == 0 .and. all(incidence == 0), &
        "periodic self-interface has zero region boundary")

    plus_region = 2
    minus_region = 1
    call initialize_region_interface_graph( &
        invalid, 1, plus_region, minus_region, status)
    call check_condition(status /= 0, &
        "out-of-range region endpoint is rejected")

    deallocate(plus_region, minus_region)
    allocate(plus_region(1), minus_region(2))
    call initialize_region_interface_graph( &
        invalid, 2, plus_region, minus_region, status)
    call check_condition(status /= 0, &
        "mismatched interface endpoint tables are rejected")

    call check_summary("region interface graph")
end program test_region_interface_graph
