module fortfem_panel_cluster_tree_3d
    !! Shared pointer-free geometric cluster tree for triangular panels.
    use fortfem_kinds, only: dp
    implicit none

    private

    public :: build_panel_cluster_tree_3d

contains

    subroutine build_panel_cluster_tree_3d( &
            vertices, triangles, leaf_size, areas, centers, permutation, &
            inverse_position, first, last, first_child, second_child, &
            node_centers, radii, node_count, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), leaf_size
        real(dp), allocatable, intent(out) :: areas(:), centers(:, :)
        integer, allocatable, intent(out) :: permutation(:)
        integer, allocatable, intent(out) :: inverse_position(:)
        integer, allocatable, intent(out) :: first(:), last(:)
        integer, allocatable, intent(out) :: first_child(:), second_child(:)
        real(dp), allocatable, intent(out) :: node_centers(:, :), radii(:)
        integer, intent(out) :: node_count, status

        integer :: node, panel, position

        status = 1
        node_count = 0
        if (size(vertices, 1) /= 3 .or. size(vertices, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (any(triangles < 1) .or. &
            any(triangles > size(vertices, 2))) return
        if (leaf_size < 1) return

        allocate( &
            areas(size(triangles, 2)), centers(3, size(triangles, 2)), &
            permutation(size(triangles, 2)), &
            inverse_position(size(triangles, 2)))
        allocate( &
            first(2*size(triangles, 2)), last(2*size(triangles, 2)), &
            first_child(2*size(triangles, 2)), &
            second_child(2*size(triangles, 2)), &
            node_centers(3, 2*size(triangles, 2)), &
            radii(2*size(triangles, 2)))
        do panel = 1, size(triangles, 2)
            call panel_geometry( &
                vertices(:, triangles(:, panel)), centers(:, panel), &
                areas(panel))
            if (areas(panel) <= 0.0_dp) return
            permutation(panel) = panel
        end do
        call build_node( &
            1, size(triangles, 2), leaf_size, vertices, triangles, centers, &
            permutation, first, last, first_child, second_child, node_centers, &
            radii, node_count, node)
        do position = 1, size(permutation)
            inverse_position(permutation(position)) = position
        end do
        status = 0
    end subroutine build_panel_cluster_tree_3d

    recursive subroutine build_node( &
            range_first, range_last, leaf_size, vertices, triangles, centers, &
            permutation, first, last, first_child, second_child, node_centers, &
            radii, node_count, node)
        integer, intent(in) :: range_first, range_last, leaf_size
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: centers(:, :)
        integer, intent(inout) :: permutation(:), first(:), last(:)
        integer, intent(inout) :: first_child(:), second_child(:), node_count
        real(dp), intent(inout) :: node_centers(:, :), radii(:)
        integer, intent(out) :: node

        real(dp) :: spans(3)
        integer :: axis, midpoint, panel, position, vertex
        real(dp) :: lower(3), upper(3)

        node_count = node_count + 1
        node = node_count
        first(node) = range_first
        last(node) = range_last
        first_child(node) = 0
        second_child(node) = 0
        lower = huge(1.0_dp)
        upper = -huge(1.0_dp)
        do position = range_first, range_last
            panel = permutation(position)
            do vertex = 1, size(triangles, 1)
                lower = min(lower, vertices(:, triangles(vertex, panel)))
                upper = max(upper, vertices(:, triangles(vertex, panel)))
            end do
        end do
        node_centers(:, node) = 0.5_dp*(lower + upper)
        radii(node) = 0.0_dp
        do position = range_first, range_last
            panel = permutation(position)
            do vertex = 1, size(triangles, 1)
                radii(node) = max(radii(node), norm2( &
                    vertices(:, triangles(vertex, panel)) - &
                    node_centers(:, node)))
            end do
        end do
        if (range_last - range_first + 1 <= leaf_size) return

        do axis = 1, 3
            spans(axis) = maxval( &
                centers(axis, permutation(range_first:range_last))) - &
                minval(centers(axis, permutation(range_first:range_last)))
        end do
        axis = maxloc(spans, dim=1)
        call sort_range( &
            range_first, range_last, axis, centers, permutation)
        midpoint = (range_first + range_last)/2
        call build_node( &
            range_first, midpoint, leaf_size, vertices, triangles, centers, &
            permutation, first, last, first_child, second_child, node_centers, &
            radii, node_count, first_child(node))
        call build_node( &
            midpoint + 1, range_last, leaf_size, vertices, triangles, centers, &
            permutation, first, last, first_child, second_child, node_centers, &
            radii, node_count, second_child(node))
    end subroutine build_node

    pure subroutine sort_range( &
            range_first, range_last, axis, centers, permutation)
        integer, intent(in) :: range_first, range_last, axis
        real(dp), intent(in) :: centers(:, :)
        integer, intent(inout) :: permutation(:)

        integer :: position, previous, selected

        do position = range_first + 1, range_last
            selected = permutation(position)
            previous = position - 1
            do while (previous >= range_first)
                if (centers(axis, permutation(previous)) <= &
                    centers(axis, selected)) exit
                permutation(previous + 1) = permutation(previous)
                previous = previous - 1
            end do
            permutation(previous + 1) = selected
        end do
    end subroutine sort_range

    pure subroutine panel_geometry(vertices, center, area)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp), intent(out) :: center(3), area

        real(dp) :: first(3), second(3)

        center = sum(vertices, dim=2)/3.0_dp
        first = vertices(:, 2) - vertices(:, 1)
        second = vertices(:, 3) - vertices(:, 1)
        area = 0.5_dp*norm2([ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)])
    end subroutine panel_geometry

end module fortfem_panel_cluster_tree_3d
