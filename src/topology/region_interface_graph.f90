module fortfem_region_interface_graph
    !! Neutral oriented graph of regions and internal interfaces.
    !!
    !! An interface column has +1 in its plus region and -1 in its minus
    !! region.  Equal endpoints are allowed for a periodic self-identification
    !! and therefore produce a zero region-boundary column.
    implicit none
    private

    type, public :: region_interface_graph_t
        private
        integer :: region_count = 0
        integer, allocatable :: plus_region(:)
        integer, allocatable :: minus_region(:)
    end type region_interface_graph_t

    public :: initialize_region_interface_graph
    public :: validate_region_interface_graph
    public :: region_interface_graph_incidence
    public :: region_interface_graph_components
    public :: region_interface_graph_cycle_basis

contains

    subroutine initialize_region_interface_graph( &
            graph, region_count, plus_region, minus_region, status)
        type(region_interface_graph_t), intent(inout) :: graph
        integer, intent(in) :: region_count
        integer, intent(in) :: plus_region(:), minus_region(:)
        integer, intent(out) :: status

        call clear_region_interface_graph(graph)
        status = 1
        if (region_count < 0) return
        if (size(plus_region) /= size(minus_region)) then
            status = 2
            return
        end if
        if (any(plus_region < 1) .or. any(plus_region > region_count) .or. &
            any(minus_region < 1) .or. any(minus_region > region_count)) then
            status = 3
            return
        end if

        graph%region_count = region_count
        allocate(graph%plus_region(size(plus_region)))
        allocate(graph%minus_region(size(minus_region)))
        graph%plus_region = plus_region
        graph%minus_region = minus_region
        status = 0
    end subroutine initialize_region_interface_graph

    subroutine validate_region_interface_graph(graph, status)
        type(region_interface_graph_t), intent(in) :: graph
        integer, intent(out) :: status

        status = 1
        if (graph%region_count < 0) return
        if (.not. allocated(graph%plus_region) .or. &
            .not. allocated(graph%minus_region)) return
        if (size(graph%plus_region) /= size(graph%minus_region)) then
            status = 2
            return
        end if
        if (any(graph%plus_region < 1) .or. &
            any(graph%plus_region > graph%region_count) .or. &
            any(graph%minus_region < 1) .or. &
            any(graph%minus_region > graph%region_count)) then
            status = 3
            return
        end if
        status = 0
    end subroutine validate_region_interface_graph

    subroutine region_interface_graph_incidence(graph, incidence, status)
        type(region_interface_graph_t), intent(in) :: graph
        integer, allocatable, intent(out) :: incidence(:, :)
        integer, intent(out) :: status
        integer :: interface

        if (allocated(incidence)) deallocate(incidence)
        call validate_region_interface_graph(graph, status)
        if (status /= 0) return
        allocate(incidence( &
            graph%region_count, size(graph%plus_region)))
        incidence = 0
        do interface = 1, size(graph%plus_region)
            incidence(graph%plus_region(interface), interface) = &
                incidence(graph%plus_region(interface), interface) + 1
            incidence(graph%minus_region(interface), interface) = &
                incidence(graph%minus_region(interface), interface) - 1
        end do
    end subroutine region_interface_graph_incidence

    subroutine region_interface_graph_components( &
            graph, components, component_count, status)
        type(region_interface_graph_t), intent(in) :: graph
        integer, allocatable, intent(out) :: components(:)
        integer, intent(out) :: component_count
        integer, intent(out) :: status
        integer :: interface, root_plus, root_minus, region

        if (allocated(components)) deallocate(components)
        component_count = 0
        call validate_region_interface_graph(graph, status)
        if (status /= 0) return
        allocate(components(graph%region_count))
        do region = 1, graph%region_count
            components(region) = region
        end do
        do interface = 1, size(graph%plus_region)
            root_plus = find_root(components, graph%plus_region(interface))
            root_minus = find_root(components, graph%minus_region(interface))
            if (root_plus /= root_minus) components(root_minus) = root_plus
        end do
        do region = 1, graph%region_count
            components(region) = find_root(components, region)
        end do
        do region = 1, graph%region_count
            if (components(region) == region) component_count = &
                component_count + 1
        end do
        call compact_component_labels(components)
        status = 0
    end subroutine region_interface_graph_components

    subroutine region_interface_graph_cycle_basis( &
            graph, cycle_basis, cycle_count, status)
        type(region_interface_graph_t), intent(in) :: graph
        integer, allocatable, intent(out) :: cycle_basis(:, :)
        integer, intent(out) :: cycle_count
        integer, intent(out) :: status
        integer, allocatable :: components(:), parent(:), parent_edge(:)
        integer, allocatable :: parent_sign(:), depth(:), queue(:)
        logical, allocatable :: visited(:), tree_edge(:)
        integer :: component_count, edge_count, region_count
        integer :: component, current, edge, head, tail, neighbor
        integer :: cycle, path_edge, first, second

        if (allocated(cycle_basis)) deallocate(cycle_basis)
        cycle_count = 0
        call validate_region_interface_graph(graph, status)
        if (status /= 0) return
        region_count = graph%region_count
        edge_count = size(graph%plus_region)
        call region_interface_graph_components( &
            graph, components, component_count, status)
        if (status /= 0) return
        cycle_count = edge_count - region_count + component_count
        if (cycle_count < 0) then
            status = 4
            cycle_count = 0
            return
        end if
        allocate(cycle_basis(edge_count, cycle_count))
        cycle_basis = 0
        if (region_count == 0) then
            status = 0
            return
        end if

        allocate(parent(region_count), parent_edge(region_count))
        allocate(parent_sign(region_count), depth(region_count))
        allocate(queue(region_count), visited(region_count))
        allocate(tree_edge(edge_count))
        visited = .false.
        tree_edge = .false.
        parent = 0
        parent_edge = 0
        parent_sign = 0
        depth = 0
        do component = 1, region_count
            if (visited(component)) cycle
            visited(component) = .true.
            parent(component) = component
            depth(component) = 0
            head = 1
            tail = 1
            queue(tail) = component
            do while (head <= tail)
                current = queue(head)
                head = head + 1
                do edge = 1, edge_count
                    if (graph%plus_region(edge) == current) then
                        neighbor = graph%minus_region(edge)
                    else if (graph%minus_region(edge) == current) then
                        neighbor = graph%plus_region(edge)
                    else
                        cycle
                    end if
                    if (visited(neighbor)) cycle
                    visited(neighbor) = .true.
                    parent(neighbor) = current
                    parent_edge(neighbor) = edge
                    if (graph%minus_region(edge) == current .and. &
                        graph%plus_region(edge) == neighbor) then
                        parent_sign(neighbor) = 1
                    else
                        parent_sign(neighbor) = -1
                    end if
                    depth(neighbor) = depth(current) + 1
                    tree_edge(edge) = .true.
                    tail = tail + 1
                    queue(tail) = neighbor
                end do
            end do
        end do

        cycle = 0
        do edge = 1, edge_count
            if (tree_edge(edge)) cycle
            cycle = cycle + 1
            cycle_basis(edge, cycle) = 1
            first = graph%plus_region(edge)
            second = graph%minus_region(edge)
            do while (depth(first) > depth(second))
                path_edge = parent_edge(first)
                cycle_basis(path_edge, cycle) = &
                    cycle_basis(path_edge, cycle) - parent_sign(first)
                first = parent(first)
            end do
            do while (depth(second) > depth(first))
                path_edge = parent_edge(second)
                cycle_basis(path_edge, cycle) = &
                    cycle_basis(path_edge, cycle) + parent_sign(second)
                second = parent(second)
            end do
            do while (first /= second)
                path_edge = parent_edge(first)
                cycle_basis(path_edge, cycle) = &
                    cycle_basis(path_edge, cycle) - parent_sign(first)
                first = parent(first)
                path_edge = parent_edge(second)
                cycle_basis(path_edge, cycle) = &
                    cycle_basis(path_edge, cycle) + parent_sign(second)
                second = parent(second)
            end do
        end do
        if (cycle /= cycle_count) then
            deallocate(cycle_basis)
            cycle_count = 0
            status = 5
            return
        end if
        status = 0
    end subroutine region_interface_graph_cycle_basis

    subroutine clear_region_interface_graph(graph)
        type(region_interface_graph_t), intent(inout) :: graph

        if (allocated(graph%plus_region)) deallocate(graph%plus_region)
        if (allocated(graph%minus_region)) deallocate(graph%minus_region)
        graph%region_count = 0
    end subroutine clear_region_interface_graph

    integer function find_root(parent, node) result(root)
        integer, intent(inout) :: parent(:)
        integer, intent(in) :: node
        integer :: current, next

        root = node
        do while (parent(root) /= root)
            root = parent(root)
        end do
        current = node
        do while (parent(current) /= current)
            next = parent(current)
            parent(current) = root
            current = next
        end do
    end function find_root

    subroutine compact_component_labels(components)
        integer, intent(inout) :: components(:)
        integer, allocatable :: labels(:)
        integer :: region, root, label

        allocate(labels(size(components)))
        labels = 0
        label = 0
        do region = 1, size(components)
            root = components(region)
            if (labels(root) == 0) then
                label = label + 1
                labels(root) = label
            end if
            components(region) = labels(root)
        end do
    end subroutine compact_component_labels

end module fortfem_region_interface_graph
