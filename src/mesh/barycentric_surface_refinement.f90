module fortfem_barycentric_surface_refinement
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: barycentric_refine_surface_mesh

contains

    subroutine barycentric_refine_surface_mesh( &
            vertices, triangles, refined_vertices, refined_triangles)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), allocatable, intent(out) :: refined_vertices(:, :)
        integer, allocatable, intent(out) :: refined_triangles(:, :)

        integer, allocatable :: edge_vertices(:, :)
        integer :: a, b, c, center, edge_count, face
        integer :: midpoint_ab, midpoint_bc, midpoint_ca

        if (size(vertices, 1) /= 3 .or. size(triangles, 1) /= 3 .or. &
            size(triangles, 2) < 1) then
            error stop "Barycentric refinement requires triangular surface"
        end if
        if (any(triangles < 1) .or. any(triangles > size(vertices, 2))) then
            error stop "Barycentric refinement has invalid vertex indices"
        end if
        allocate(edge_vertices(2, 3*size(triangles, 2)))
        edge_count = 0
        do face = 1, size(triangles, 2)
            call register_edge( &
                triangles(1, face), triangles(2, face), edge_vertices, &
                edge_count)
            call register_edge( &
                triangles(2, face), triangles(3, face), edge_vertices, &
                edge_count)
            call register_edge( &
                triangles(3, face), triangles(1, face), edge_vertices, &
                edge_count)
        end do
        allocate(refined_vertices(3, &
            size(vertices, 2) + edge_count + size(triangles, 2)))
        refined_vertices(:, :size(vertices, 2)) = vertices
        do a = 1, edge_count
            refined_vertices(:, size(vertices, 2) + a) = 0.5_dp*( &
                vertices(:, edge_vertices(1, a)) + &
                vertices(:, edge_vertices(2, a)))
        end do
        do face = 1, size(triangles, 2)
            refined_vertices(:, size(vertices, 2) + edge_count + face) = &
                sum(vertices(:, triangles(:, face)), dim=2)/3.0_dp
        end do
        allocate(refined_triangles(3, 6*size(triangles, 2)))
        do face = 1, size(triangles, 2)
            a = triangles(1, face)
            b = triangles(2, face)
            c = triangles(3, face)
            midpoint_ab = size(vertices, 2) + &
                find_edge(a, b, edge_vertices, edge_count)
            midpoint_bc = size(vertices, 2) + &
                find_edge(b, c, edge_vertices, edge_count)
            midpoint_ca = size(vertices, 2) + &
                find_edge(c, a, edge_vertices, edge_count)
            center = size(vertices, 2) + edge_count + face
            refined_triangles(:, 6*face - 5) = [a, midpoint_ab, center]
            refined_triangles(:, 6*face - 4) = [a, center, midpoint_ca]
            refined_triangles(:, 6*face - 3) = [b, midpoint_bc, center]
            refined_triangles(:, 6*face - 2) = [b, center, midpoint_ab]
            refined_triangles(:, 6*face - 1) = [c, midpoint_ca, center]
            refined_triangles(:, 6*face) = [c, center, midpoint_bc]
        end do
    end subroutine barycentric_refine_surface_mesh

    subroutine register_edge(first, second, edges, edge_count)
        integer, intent(in) :: first, second
        integer, intent(inout) :: edges(:, :), edge_count

        integer :: edge, ordered(2)

        ordered = [min(first, second), max(first, second)]
        do edge = 1, edge_count
            if (all(edges(:, edge) == ordered)) return
        end do
        edge_count = edge_count + 1
        edges(:, edge_count) = ordered
    end subroutine register_edge

    pure integer function find_edge( &
            first, second, edges, edge_count) result(index)
        integer, intent(in) :: first, second, edges(:, :), edge_count

        integer :: edge, ordered(2)

        ordered = [min(first, second), max(first, second)]
        index = 0
        do edge = 1, edge_count
            if (all(edges(:, edge) == ordered)) then
                index = edge
                return
            end if
        end do
    end function find_edge

end module fortfem_barycentric_surface_refinement
