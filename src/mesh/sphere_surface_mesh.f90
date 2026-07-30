module fortfem_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: generate_sphere_surface_mesh

contains

    subroutine generate_sphere_surface_mesh( &
            radius, refinement_level, vertices, triangles)
        real(dp), intent(in) :: radius
        integer, intent(in) :: refinement_level
        real(dp), allocatable, intent(out) :: vertices(:, :)
        integer, allocatable, intent(out) :: triangles(:, :)

        integer :: refinement

        if (radius <= 0.0_dp .or. refinement_level < 0) then
            error stop "Sphere mesh requires positive radius and refinement"
        end if
        allocate(vertices(3, 6), triangles(3, 8))
        vertices(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
        vertices(:, 2) = [-1.0_dp, 0.0_dp, 0.0_dp]
        vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        vertices(:, 4) = [0.0_dp, -1.0_dp, 0.0_dp]
        vertices(:, 5) = [0.0_dp, 0.0_dp, 1.0_dp]
        vertices(:, 6) = [0.0_dp, 0.0_dp, -1.0_dp]
        triangles(:, 1) = [1, 3, 5]
        triangles(:, 2) = [3, 2, 5]
        triangles(:, 3) = [2, 4, 5]
        triangles(:, 4) = [4, 1, 5]
        triangles(:, 5) = [3, 1, 6]
        triangles(:, 6) = [2, 3, 6]
        triangles(:, 7) = [4, 2, 6]
        triangles(:, 8) = [1, 4, 6]
        do refinement = 1, refinement_level
            call refine_unit_sphere(vertices, triangles)
        end do
        vertices = radius*vertices
    end subroutine generate_sphere_surface_mesh

    subroutine refine_unit_sphere(vertices, triangles)
        real(dp), allocatable, intent(inout) :: vertices(:, :)
        integer, allocatable, intent(inout) :: triangles(:, :)

        real(dp), allocatable :: expanded_vertices(:, :)
        integer, allocatable :: edge_midpoints(:), edge_vertices(:, :)
        integer, allocatable :: refined(:, :)
        integer :: edge_count, midpoint(3), old_vertex_count
        integer :: triangle, vertex_count

        old_vertex_count = size(vertices, 2)
        allocate(expanded_vertices(3, old_vertex_count + 3*size(triangles, 2)))
        expanded_vertices(:, :old_vertex_count) = vertices
        allocate(edge_vertices(2, 3*size(triangles, 2)))
        allocate(edge_midpoints(3*size(triangles, 2)))
        allocate(refined(3, 4*size(triangles, 2)))
        vertex_count = old_vertex_count
        edge_count = 0
        do triangle = 1, size(triangles, 2)
            midpoint(1) = midpoint_vertex( &
                triangles(1, triangle), triangles(2, triangle), vertices, &
                expanded_vertices, edge_vertices, edge_midpoints, &
                edge_count, vertex_count)
            midpoint(2) = midpoint_vertex( &
                triangles(2, triangle), triangles(3, triangle), vertices, &
                expanded_vertices, edge_vertices, edge_midpoints, &
                edge_count, vertex_count)
            midpoint(3) = midpoint_vertex( &
                triangles(3, triangle), triangles(1, triangle), vertices, &
                expanded_vertices, edge_vertices, edge_midpoints, &
                edge_count, vertex_count)
            refined(:, 4*triangle - 3) = &
                [triangles(1, triangle), midpoint(1), midpoint(3)]
            refined(:, 4*triangle - 2) = &
                [midpoint(1), triangles(2, triangle), midpoint(2)]
            refined(:, 4*triangle - 1) = &
                [midpoint(3), midpoint(2), triangles(3, triangle)]
            refined(:, 4*triangle) = midpoint
        end do
        vertices = expanded_vertices(:, :vertex_count)
        call move_alloc(refined, triangles)
    end subroutine refine_unit_sphere

    function midpoint_vertex( &
            first, second, vertices, expanded_vertices, edge_vertices, &
            edge_midpoints, edge_count, vertex_count) result(midpoint)
        integer, intent(in) :: first, second
        real(dp), intent(in) :: vertices(:, :)
        real(dp), intent(inout) :: expanded_vertices(:, :)
        integer, intent(inout) :: edge_vertices(:, :), edge_midpoints(:)
        integer, intent(inout) :: edge_count, vertex_count
        integer :: midpoint

        integer :: edge, ordered(2)
        real(dp) :: midpoint_coordinate(3)

        ordered = [min(first, second), max(first, second)]
        do edge = 1, edge_count
            if (all(edge_vertices(:, edge) == ordered)) then
                midpoint = edge_midpoints(edge)
                return
            end if
        end do
        edge_count = edge_count + 1
        edge_vertices(:, edge_count) = ordered
        vertex_count = vertex_count + 1
        midpoint_coordinate = vertices(:, first) + vertices(:, second)
        expanded_vertices(:, vertex_count) = &
            midpoint_coordinate/norm2(midpoint_coordinate)
        edge_midpoints(edge_count) = vertex_count
        midpoint = vertex_count
    end function midpoint_vertex

end module fortfem_sphere_surface_mesh
