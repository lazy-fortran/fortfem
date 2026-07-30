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
    use fortfem_panel_cluster_tree_3d, only: build_panel_cluster_tree_3d
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
        integer :: local_triangle(3, 1), node, node_count, panel

        status = 1
        interaction_count = 0
        if (size(vertices, 1) /= 3 .or. size(vertices, 2) < 3) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (size(density) /= size(triangles, 2)) return
        if (any(triangles < 1) .or. &
            any(triangles > size(vertices, 2))) return
        if (opening_angle <= 0.0_dp .or. opening_angle >= 1.0_dp) return
        if (leaf_size < 1) return

        call build_panel_cluster_tree_3d( &
            vertices, triangles, leaf_size, areas, centers, permutation, &
            inverse_position, first, last, first_child, second_child, &
            node_centers, radii, node_count, status)
        if (status /= 0) return
        allocate(result(size(density)), charges(2*size(density)))
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

end module fortfem_laplace_hierarchical_3d
