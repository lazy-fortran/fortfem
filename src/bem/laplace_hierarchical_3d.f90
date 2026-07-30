module fortfem_laplace_hierarchical_3d
    !! Pointer-free Barnes-Hut application of the P0 Laplace single layer.
    !!
    !! Analytical Galerkin self panels are retained. Far clusters use their
    !! area-weighted monopole; near leaves use panel-centroid interactions.
    !! Storage is O(N), and accepted interaction counts are O(N log N) for
    !! quasi-uniform closed surfaces at fixed opening angle.
    use fortfem_kinds, only: dp
    use fortfem_laplace_galerkin_3d, only: &
        assemble_laplace_single_layer_p0_3d
    implicit none

    private

    public :: apply_laplace_single_layer_p0_hierarchical_3d

contains

    subroutine apply_laplace_single_layer_p0_hierarchical_3d( &
            vertices, triangles, density, opening_angle, leaf_size, result, &
            status, interaction_count)
        real(dp), intent(in) :: vertices(:, :), density(:), opening_angle
        integer, intent(in) :: triangles(:, :), leaf_size
        real(dp), allocatable, intent(out) :: result(:)
        integer, intent(out) :: status, interaction_count

        integer, allocatable :: first_child(:), inverse_position(:)
        integer, allocatable :: last(:), first(:), permutation(:)
        integer, allocatable :: second_child(:)
        real(dp), allocatable :: areas(:), centers(:, :), charges(:)
        real(dp), allocatable :: node_centers(:, :), radii(:)
        real(dp), allocatable :: self_matrix(:, :)
        real(dp) :: local_vertices(3, 3)
        integer :: local_triangle(3, 1), node, node_count, panel, position

        status = 1
        interaction_count = 0
        if (size(vertices, 1) /= 3 .or. size(vertices, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (size(density) /= size(triangles, 2)) return
        if (any(triangles < 1) .or. &
            any(triangles > size(vertices, 2))) return
        if (opening_angle <= 0.0_dp .or. opening_angle >= 1.0_dp) return
        if (leaf_size < 1) return

        allocate( &
            result(size(density)), areas(size(density)), &
            centers(3, size(density)), permutation(size(density)), &
            inverse_position(size(density)))
        allocate( &
            first(2*size(density)), last(2*size(density)), &
            first_child(2*size(density)), second_child(2*size(density)), &
            node_centers(3, 2*size(density)), radii(2*size(density)), &
            charges(2*size(density)))
        do panel = 1, size(density)
            call panel_geometry( &
                vertices(:, triangles(:, panel)), centers(:, panel), &
                areas(panel))
            if (areas(panel) <= 0.0_dp) return
            permutation(panel) = panel
        end do
        node_count = 0
        call build_node( &
            1, size(density), leaf_size, centers, permutation, first, last, &
            first_child, second_child, node_centers, radii, node_count, node)
        do position = 1, size(permutation)
            inverse_position(permutation(position)) = position
        end do
        do node = node_count, 1, -1
            charges(node) = sum( &
                density(permutation(first(node):last(node)))* &
                areas(permutation(first(node):last(node))))
        end do

        result = 0.0_dp
        local_triangle(:, 1) = [1, 2, 3]
        do panel = 1, size(density)
            local_vertices = vertices(:, triangles(:, panel))
            call assemble_laplace_single_layer_p0_3d( &
                local_vertices, local_triangle, 8, self_matrix, status)
            if (status /= 0) return
            result(panel) = self_matrix(1, 1)*density(panel)
            interaction_count = interaction_count + 1
            call apply_node( &
                1, panel, inverse_position(panel), opening_angle, &
                centers, areas, density, permutation, first, last, &
                first_child, second_child, node_centers, radii, charges, &
                result(panel), interaction_count)
        end do
        status = 0
    end subroutine apply_laplace_single_layer_p0_hierarchical_3d

    recursive subroutine build_node( &
            range_first, range_last, leaf_size, centers, permutation, first, &
            last, first_child, second_child, node_centers, radii, node_count, &
            node)
        integer, intent(in) :: range_first, range_last, leaf_size
        real(dp), intent(in) :: centers(:, :)
        integer, intent(inout) :: permutation(:), first(:), last(:)
        integer, intent(inout) :: first_child(:), second_child(:), node_count
        real(dp), intent(inout) :: node_centers(:, :), radii(:)
        integer, intent(out) :: node

        real(dp) :: spans(3)
        integer :: axis, midpoint, position

        node_count = node_count + 1
        node = node_count
        first(node) = range_first
        last(node) = range_last
        first_child(node) = 0
        second_child(node) = 0
        node_centers(:, node) = 0.0_dp
        do position = range_first, range_last
            node_centers(:, node) = node_centers(:, node) + &
                centers(:, permutation(position))
        end do
        node_centers(:, node) = node_centers(:, node)/ &
            real(range_last - range_first + 1, dp)
        radii(node) = 0.0_dp
        do position = range_first, range_last
            radii(node) = max(radii(node), norm2( &
                centers(:, permutation(position)) - node_centers(:, node)))
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
            range_first, midpoint, leaf_size, centers, permutation, first, &
            last, first_child, second_child, node_centers, radii, node_count, &
            first_child(node))
        call build_node( &
            midpoint + 1, range_last, leaf_size, centers, permutation, first, &
            last, first_child, second_child, node_centers, radii, node_count, &
            second_child(node))
    end subroutine build_node

    recursive subroutine apply_node( &
            node, target, target_position, opening_angle, centers, areas, &
            density, permutation, first, last, first_child, second_child, &
            node_centers, radii, charges, value, interaction_count)
        integer, intent(in) :: node, target, target_position
        real(dp), intent(in) :: opening_angle, centers(:, :), areas(:)
        real(dp), intent(in) :: density(:), node_centers(:, :), radii(:)
        real(dp), intent(in) :: charges(:)
        integer, intent(in) :: permutation(:), first(:), last(:)
        integer, intent(in) :: first_child(:), second_child(:)
        real(dp), intent(inout) :: value
        integer, intent(inout) :: interaction_count

        real(dp) :: distance
        integer :: panel, position
        logical :: contains_target

        contains_target = target_position >= first(node) .and. &
            target_position <= last(node)
        distance = norm2(centers(:, target) - node_centers(:, node))
        if (.not. contains_target .and. distance > 0.0_dp .and. &
            radii(node)/distance < opening_angle) then
            value = value + areas(target)*charges(node)/ &
                (4.0_dp*acos(-1.0_dp)*distance)
            interaction_count = interaction_count + 1
            return
        end if
        if (first_child(node) == 0) then
            do position = first(node), last(node)
                panel = permutation(position)
                if (panel == target) cycle
                distance = norm2(centers(:, target) - centers(:, panel))
                value = value + areas(target)*areas(panel)*density(panel)/ &
                    (4.0_dp*acos(-1.0_dp)*distance)
                interaction_count = interaction_count + 1
            end do
            return
        end if
        call apply_node( &
            first_child(node), target, target_position, opening_angle, &
            centers, areas, density, permutation, first, last, first_child, &
            second_child, node_centers, radii, charges, value, &
            interaction_count)
        call apply_node( &
            second_child(node), target, target_position, opening_angle, &
            centers, areas, density, permutation, first, last, first_child, &
            second_child, node_centers, radii, charges, value, &
            interaction_count)
    end subroutine apply_node

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

end module fortfem_laplace_hierarchical_3d
