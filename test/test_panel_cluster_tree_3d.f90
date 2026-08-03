program test_panel_cluster_tree_3d
    !! Independent geometry oracle for the internal BEM cluster tree.
    !!
    !! A node radius must enclose every panel vertex, not only panel
    !! centroids.  Otherwise an admissible far interaction can skip part of
    !! a large triangle.  The check below is deliberately independent of the
    !! tree implementation: it walks every node range and recomputes the
    !! vertex distances from the returned node centre.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_panel_cluster_tree_3d, only: build_panel_cluster_tree_3d
    implicit none

    real(dp) :: vertices(3, 12)
    real(dp), allocatable :: areas(:), centers(:, :), radii(:)
    integer :: triangles(3, 4)
    integer, allocatable :: permutation(:), inverse_position(:)
    integer, allocatable :: first(:), last(:), first_child(:), second_child(:)
    real(dp), allocatable :: node_centers(:, :)
    real(dp) :: expected_area, distance
    integer :: node_count, status, node, position, panel, vertex
    logical :: all_passed

    ! Four nondegenerate panels with deliberately large vertex offsets from
    ! their centroids.  Their supports must be enclosed by every ancestor.
    vertices = reshape([ &
        -2.0_dp, -1.0_dp,  0.0_dp,  2.0_dp, -1.0_dp,  0.0_dp, &
         0.0_dp,  3.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  4.0_dp, &
         8.0_dp, -1.0_dp,  0.0_dp, 12.0_dp, -1.0_dp,  0.0_dp, &
        10.0_dp,  3.0_dp,  0.0_dp, 10.0_dp,  0.0_dp,  4.0_dp, &
        10.0_dp,  0.0_dp,  0.0_dp, 10.0_dp,  4.0_dp,  4.0_dp, &
        10.0_dp,  0.0_dp,  4.0_dp, 10.0_dp,  0.0_dp,  0.0_dp], [3, 12])
    triangles = reshape([ &
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [3, 4])

    all_passed = .true.
    call build_panel_cluster_tree_3d( &
        vertices, triangles, 2, areas, centers, permutation, inverse_position, &
        first, last, first_child, second_child, node_centers, radii, &
        node_count, status)
    call record_condition(status == 0, "cluster tree accepts valid panels")
    call record_condition(node_count >= 3, "cluster tree creates root and leaves")
    call record_condition(size(areas) == 4 .and. size(centers, 2) == 4, &
        "cluster tree returns one geometry record per panel")

    do panel = 1, size(areas)
        expected_area = 0.5_dp*sqrt(sum(cross_product( &
            vertices(:, triangles(2, panel)) - vertices(:, triangles(1, panel)), &
            vertices(:, triangles(3, panel)) - vertices(:, triangles(1, panel)))**2))
        call record_condition(abs(areas(panel) - expected_area) < 1.0e-13_dp, &
            "panel area matches independent cross-product oracle")
    end do
    do panel = 1, size(permutation)
        call record_condition(inverse_position(permutation(panel)) == panel, &
            "permutation and inverse permutation are consistent")
    end do

    do node = 1, node_count
        do position = first(node), last(node)
            panel = permutation(position)
            do vertex = 1, 3
                distance = sqrt(sum(( &
                    vertices(:, triangles(vertex, panel)) - node_centers(:, node))**2))
                call record_condition(distance <= radii(node) + 1.0e-13_dp, &
                    "node radius encloses every panel vertex")
            end do
        end do
    end do

    call build_panel_cluster_tree_3d( &
        vertices, triangles, 0, areas, centers, permutation, inverse_position, &
        first, last, first_child, second_child, node_centers, radii, &
        node_count, status)
    call record_condition(status /= 0, "cluster tree rejects zero leaf size")
    call check_summary("3-D BEM panel cluster tree")
    if (.not. all_passed) error stop 1

contains

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_panel_cluster_tree_3d
