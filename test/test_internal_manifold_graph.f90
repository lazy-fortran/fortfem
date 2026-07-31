program test_internal_manifold_graph
    use check, only: check_condition, check_summary
    use fortfem_api, only: initialize_internal_manifold_graph, &
        internal_manifold_graph_closed, internal_manifold_graph_components, &
        internal_manifold_graph_junction_incidence, &
        internal_manifold_graph_region_incidence, internal_manifold_graph_t, &
        validate_internal_manifold_graph
    implicit none

    type(internal_manifold_graph_t) :: slab, cylinder, sphere, torus, invalid
    integer, allocatable :: region_incidence(:, :), junction_incidence(:, :)
    integer, allocatable :: components(:)
    logical, allocatable :: closed(:)
    integer :: component_count, status

    call initialize_internal_manifold_graph( &
        slab, 3, 3, [2, 3], [1, 2], [1, 2], [2, 3], status)
    call check_condition(status == 0, "open slab manifold chain initializes")
    call validate_internal_manifold_graph(slab, status)
    call check_condition(status == 0, "open slab manifold chain validates")
    call internal_manifold_graph_region_incidence( &
        slab, region_incidence, status)
    call check_condition(status == 0 .and. all(region_incidence == reshape([ &
        -1, 1, 0, 0, -1, 1], [3, 2])), &
        "manifold region incidence preserves plus-minus orientation")
    call internal_manifold_graph_junction_incidence( &
        slab, junction_incidence, status)
    call check_condition(status == 0 .and. all(junction_incidence == reshape([ &
        -1, 1, 0, 0, -1, 1], [3, 2])), &
        "open manifold junction incidence preserves boundary orientation")
    call internal_manifold_graph_closed(slab, closed, status)
    call check_condition(status == 0 .and. all(.not. closed), &
        "slab interfaces with distinct endpoints are open")
    call internal_manifold_graph_components( &
        slab, components, component_count, status)
    call check_condition(status == 0 .and. component_count == 1 .and. &
        all(components == [1, 1]), &
        "connected slab manifold chain has one component")

    call initialize_internal_manifold_graph( &
        cylinder, 1, 1, [1], [1], [1], [1], status)
    call internal_manifold_graph_junction_incidence( &
        cylinder, junction_incidence, status)
    call internal_manifold_graph_closed(cylinder, closed, status)
    call check_condition(status == 0 .and. all(junction_incidence == 0) .and. &
        closed(1), "periodic cylinder self-identification is closed")

    call initialize_internal_manifold_graph( &
        sphere, 1, 0, [1], [1], [0], [0], status)
    call internal_manifold_graph_junction_incidence( &
        sphere, junction_incidence, status)
    call internal_manifold_graph_closed(sphere, closed, status)
    call check_condition(status == 0 .and. size(junction_incidence, 1) == 0 &
        .and. closed(1), "boundaryless sphere manifold has no junction flux")

    call initialize_internal_manifold_graph( &
        torus, 2, 1, [2, 2], [1, 1], [1, 1], [1, 1], status)
    call internal_manifold_graph_components( &
        torus, components, component_count, status)
    call internal_manifold_graph_closed(torus, closed, status)
    call check_condition(status == 0 .and. component_count == 1 .and. &
        all(components == [1, 1]) .and. all(closed), &
        "toroidal periodic manifolds share one closed component")

    call initialize_internal_manifold_graph( &
        invalid, 1, 1, [1], [1], [0], [1], status)
    call check_condition(status /= 0, &
        "a manifold cannot have only one boundary endpoint")
    call initialize_internal_manifold_graph( &
        invalid, 1, 1, [1], [1], [1], [2], status)
    call check_condition(status /= 0, "out-of-range junction is rejected")

    call check_summary("internal manifold graph")
end program test_internal_manifold_graph
