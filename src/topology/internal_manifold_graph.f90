module fortfem_internal_manifold_graph
    !! Neutral topology for oriented internal manifolds.
    !!
    !! A manifold has a plus and minus region and optional oriented boundary
    !! junctions.  Endpoint zero denotes a boundaryless closed manifold;
    !! equal nonzero endpoints denote a periodic self-identification.  This
    !! module stores no geometry, trace space, constitutive law, or mesh.
    implicit none
    private

    type, public :: internal_manifold_graph_t
        private
        integer :: region_count = 0
        integer :: manifold_count = 0
        integer :: junction_count = 0
        integer, allocatable :: plus_region(:)
        integer, allocatable :: minus_region(:)
        integer, allocatable :: start_junction(:)
        integer, allocatable :: end_junction(:)
    end type internal_manifold_graph_t

    public :: initialize_internal_manifold_graph
    public :: validate_internal_manifold_graph
    public :: internal_manifold_graph_region_incidence
    public :: internal_manifold_graph_junction_incidence
    public :: internal_manifold_graph_closed
    public :: internal_manifold_graph_components

contains

    subroutine initialize_internal_manifold_graph( &
            graph, region_count, junction_count, plus_region, minus_region, &
            start_junction, end_junction, status)
        type(internal_manifold_graph_t), intent(inout) :: graph
        integer, intent(in) :: region_count, junction_count
        integer, intent(in) :: plus_region(:), minus_region(:)
        integer, intent(in) :: start_junction(:), end_junction(:)
        integer, intent(out) :: status
        integer :: manifold_count

        call clear_internal_manifold_graph(graph)
        status = 1
        manifold_count = size(plus_region)
        if (region_count < 0 .or. junction_count < 0) return
        if (size(minus_region) /= manifold_count .or. &
            size(start_junction) /= manifold_count .or. &
            size(end_junction) /= manifold_count) return
        if (any(plus_region < 1) .or. any(plus_region > region_count) .or. &
            any(minus_region < 1) .or. any(minus_region > region_count)) return
        if (any(start_junction < 0) .or. &
            any(start_junction > junction_count) .or. &
            any(end_junction < 0) .or. any(end_junction > junction_count)) return
        if (any((start_junction == 0) .neqv. (end_junction == 0))) return

        graph%region_count = region_count
        graph%manifold_count = manifold_count
        graph%junction_count = junction_count
        allocate(graph%plus_region(manifold_count))
        allocate(graph%minus_region(manifold_count))
        allocate(graph%start_junction(manifold_count))
        allocate(graph%end_junction(manifold_count))
        graph%plus_region = plus_region
        graph%minus_region = minus_region
        graph%start_junction = start_junction
        graph%end_junction = end_junction
        status = 0
    end subroutine initialize_internal_manifold_graph

    subroutine validate_internal_manifold_graph(graph, status)
        type(internal_manifold_graph_t), intent(in) :: graph
        integer, intent(out) :: status
        integer :: manifold_count

        status = 1
        if (graph%region_count < 0 .or. graph%manifold_count < 0 .or. &
            graph%junction_count < 0) return
        if (.not. allocated(graph%plus_region) .or. &
            .not. allocated(graph%minus_region) .or. &
            .not. allocated(graph%start_junction) .or. &
            .not. allocated(graph%end_junction)) return
        manifold_count = graph%manifold_count
        if (size(graph%plus_region) /= manifold_count .or. &
            size(graph%minus_region) /= manifold_count .or. &
            size(graph%start_junction) /= manifold_count .or. &
            size(graph%end_junction) /= manifold_count) return
        if (any(graph%plus_region < 1) .or. &
            any(graph%plus_region > graph%region_count) .or. &
            any(graph%minus_region < 1) .or. &
            any(graph%minus_region > graph%region_count)) return
        if (any(graph%start_junction < 0) .or. &
            any(graph%start_junction > graph%junction_count) .or. &
            any(graph%end_junction < 0) .or. &
            any(graph%end_junction > graph%junction_count)) return
        if (any((graph%start_junction == 0) .neqv. &
            (graph%end_junction == 0))) return
        status = 0
    end subroutine validate_internal_manifold_graph

    subroutine internal_manifold_graph_region_incidence( &
            graph, incidence, status)
        type(internal_manifold_graph_t), intent(in) :: graph
        integer, allocatable, intent(out) :: incidence(:, :)
        integer, intent(out) :: status
        integer :: manifold

        if (allocated(incidence)) deallocate(incidence)
        call validate_internal_manifold_graph(graph, status)
        if (status /= 0) return
        allocate(incidence(graph%region_count, graph%manifold_count))
        incidence = 0
        do manifold = 1, graph%manifold_count
            incidence(graph%plus_region(manifold), manifold) = &
                incidence(graph%plus_region(manifold), manifold) + 1
            incidence(graph%minus_region(manifold), manifold) = &
                incidence(graph%minus_region(manifold), manifold) - 1
        end do
        status = 0
    end subroutine internal_manifold_graph_region_incidence

    subroutine internal_manifold_graph_junction_incidence( &
            graph, incidence, status)
        type(internal_manifold_graph_t), intent(in) :: graph
        integer, allocatable, intent(out) :: incidence(:, :)
        integer, intent(out) :: status
        integer :: manifold, start, finish

        if (allocated(incidence)) deallocate(incidence)
        call validate_internal_manifold_graph(graph, status)
        if (status /= 0) return
        allocate(incidence(graph%junction_count, graph%manifold_count))
        incidence = 0
        do manifold = 1, graph%manifold_count
            start = graph%start_junction(manifold)
            finish = graph%end_junction(manifold)
            if (start > 0) incidence(start, manifold) = &
                incidence(start, manifold) - 1
            if (finish > 0) incidence(finish, manifold) = &
                incidence(finish, manifold) + 1
        end do
        status = 0
    end subroutine internal_manifold_graph_junction_incidence

    subroutine internal_manifold_graph_closed(graph, closed, status)
        type(internal_manifold_graph_t), intent(in) :: graph
        logical, allocatable, intent(out) :: closed(:)
        integer, intent(out) :: status

        if (allocated(closed)) deallocate(closed)
        call validate_internal_manifold_graph(graph, status)
        if (status /= 0) then
            allocate(closed(0))
            return
        end if
        allocate(closed(graph%manifold_count))
        closed = graph%start_junction == graph%end_junction
        status = 0
    end subroutine internal_manifold_graph_closed

    subroutine internal_manifold_graph_components( &
            graph, components, component_count, status)
        type(internal_manifold_graph_t), intent(in) :: graph
        integer, allocatable, intent(out) :: components(:)
        integer, intent(out) :: component_count
        integer, intent(out) :: status
        integer, allocatable :: parent(:)
        integer :: manifold, other, root_first, root_second

        if (allocated(components)) deallocate(components)
        component_count = 0
        call validate_internal_manifold_graph(graph, status)
        if (status /= 0) return
        allocate(components(graph%manifold_count))
        allocate(parent(graph%manifold_count))
        do manifold = 1, graph%manifold_count
            parent(manifold) = manifold
        end do
        do manifold = 1, graph%manifold_count
            do other = manifold + 1, graph%manifold_count
                if (.not. share_junction(graph, manifold, other)) cycle
                root_first = find_root(parent, manifold)
                root_second = find_root(parent, other)
                if (root_first /= root_second) parent(root_second) = root_first
            end do
        end do
        do manifold = 1, graph%manifold_count
            components(manifold) = find_root(parent, manifold)
        end do
        call compact_component_labels(components, component_count)
        status = 0
    end subroutine internal_manifold_graph_components

    logical function share_junction(graph, first, second) result(shared)
        type(internal_manifold_graph_t), intent(in) :: graph
        integer, intent(in) :: first, second

        shared = .false.
        if (graph%start_junction(first) > 0) then
            shared = graph%start_junction(first) == graph%start_junction(second) &
                .or. graph%start_junction(first) == graph%end_junction(second)
        end if
        if (.not. shared .and. graph%end_junction(first) > 0) then
            shared = graph%end_junction(first) == graph%start_junction(second) &
                .or. graph%end_junction(first) == graph%end_junction(second)
        end if
    end function share_junction

    subroutine clear_internal_manifold_graph(graph)
        type(internal_manifold_graph_t), intent(inout) :: graph

        if (allocated(graph%plus_region)) deallocate(graph%plus_region)
        if (allocated(graph%minus_region)) deallocate(graph%minus_region)
        if (allocated(graph%start_junction)) deallocate(graph%start_junction)
        if (allocated(graph%end_junction)) deallocate(graph%end_junction)
        graph%region_count = 0
        graph%manifold_count = 0
        graph%junction_count = 0
    end subroutine clear_internal_manifold_graph

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

    subroutine compact_component_labels(components, component_count)
        integer, intent(inout) :: components(:)
        integer, intent(out) :: component_count
        integer, allocatable :: labels(:)
        integer :: item, root

        allocate(labels(size(components)))
        labels = 0
        component_count = 0
        do item = 1, size(components)
            root = components(item)
            if (labels(root) == 0) then
                component_count = component_count + 1
                labels(root) = component_count
            end if
            components(item) = labels(root)
        end do
    end subroutine compact_component_labels

end module fortfem_internal_manifold_graph
