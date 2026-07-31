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
