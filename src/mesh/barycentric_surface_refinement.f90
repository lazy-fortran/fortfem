module fortfem_barycentric_surface_refinement
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: barycentric_refine_surface_mesh
    public :: barycentric_refine_torus_surface_mesh

contains

    subroutine barycentric_refine_torus_surface_mesh( &
            vertices, triangles, parameters, major_radius, minor_radius, &
            refined_vertices, refined_triangles, refined_parameters, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: triangles(:, :)
        real(dp), allocatable, intent(out) :: refined_vertices(:, :)
        integer, allocatable, intent(out) :: refined_triangles(:, :)
        real(dp), allocatable, intent(out) :: refined_parameters(:, :)
        integer, intent(out) :: status

        integer, allocatable :: centroid_vertices(:), edge_midpoints(:)
        integer, allocatable :: primal_edges(:, :)
        real(dp) :: local_parameters(2, 3), theta, phi
        integer :: edge, face, local, vertex

        status = 1
        if (allocated(refined_vertices)) deallocate(refined_vertices)
        if (allocated(refined_triangles)) deallocate(refined_triangles)
        if (allocated(refined_parameters)) deallocate(refined_parameters)
        if (size(vertices, 1) /= 3 .or. size(parameters, 1) /= 2) return
        if (size(parameters, 2) /= size(vertices, 2)) return
        if (size(triangles, 1) /= 3 .or. size(triangles, 2) < 1) return
        if (any(triangles < 1) .or. any(triangles > size(vertices, 2))) return
        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) return
        call barycentric_refine_surface_mesh( &
            vertices, triangles, refined_vertices, refined_triangles, &
            primal_edges, edge_midpoints, centroid_vertices)
        allocate(refined_parameters(2, size(refined_vertices, 2)))
        refined_parameters(:, :size(parameters, 2)) = parameters
        do edge = 1, size(primal_edges, 2)
            local_parameters(:, 1) = parameters(:, primal_edges(1, edge))
            local_parameters(:, 2) = parameters(:, primal_edges(2, edge))
            call unwrap_parameters( &
                local_parameters(:, 1), local_parameters(:, 2))
            refined_parameters(:, edge_midpoints(edge)) = &
                0.5_dp*(local_parameters(:, 1) + local_parameters(:, 2))
        end do
        do face = 1, size(triangles, 2)
            do local = 1, 3
                local_parameters(:, local) = &
                    parameters(:, triangles(local, face))
            end do
            call unwrap_parameters( &
                local_parameters(:, 1), local_parameters(:, 2))
            call unwrap_parameters( &
                local_parameters(:, 1), local_parameters(:, 3))
            refined_parameters(:, centroid_vertices(face)) = &
                sum(local_parameters, dim=2)/3.0_dp
        end do
        do vertex = 1, size(refined_vertices, 2)
            theta = refined_parameters(1, vertex)
            phi = refined_parameters(2, vertex)
            refined_vertices(:, vertex) = [ &
                (major_radius + minor_radius*cos(theta))*cos(phi), &
                (major_radius + minor_radius*cos(theta))*sin(phi), &
                minor_radius*sin(theta)]
        end do
        status = 0
    end subroutine barycentric_refine_torus_surface_mesh

    pure subroutine unwrap_parameters(reference, value)
        real(dp), intent(in) :: reference(2)
        real(dp), intent(inout) :: value(2)

        integer :: coordinate

        do coordinate = 1, 2
            do while (value(coordinate) - reference(coordinate) > &
                    acos(-1.0_dp))
                value(coordinate) = value(coordinate) - 2.0_dp*acos(-1.0_dp)
            end do
            do while (value(coordinate) - reference(coordinate) < &
                    -acos(-1.0_dp))
                value(coordinate) = value(coordinate) + 2.0_dp*acos(-1.0_dp)
            end do
        end do
    end subroutine unwrap_parameters

    subroutine barycentric_refine_surface_mesh( &
            vertices, triangles, refined_vertices, refined_triangles, &
            primal_edges, edge_midpoint_vertices, face_centroid_vertices, &
            parent_faces, local_sectors)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), allocatable, intent(out) :: refined_vertices(:, :)
        integer, allocatable, intent(out) :: refined_triangles(:, :)
        integer, allocatable, intent(out), optional :: primal_edges(:, :)
        integer, allocatable, intent(out), optional :: edge_midpoint_vertices(:)
        integer, allocatable, intent(out), optional :: face_centroid_vertices(:)
        integer, allocatable, intent(out), optional :: parent_faces(:)
        integer, allocatable, intent(out), optional :: local_sectors(:)

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
        if (present(primal_edges)) then
            allocate(primal_edges(2, edge_count))
            primal_edges = edge_vertices(:, :edge_count)
        end if
        if (present(edge_midpoint_vertices)) then
            allocate(edge_midpoint_vertices(edge_count))
            edge_midpoint_vertices = [ &
                (size(vertices, 2) + a, a=1, edge_count)]
        end if
        if (present(face_centroid_vertices)) then
            allocate(face_centroid_vertices(size(triangles, 2)))
            face_centroid_vertices = [ &
                (size(vertices, 2) + edge_count + face, &
                face=1, size(triangles, 2))]
        end if
        allocate(refined_triangles(3, 6*size(triangles, 2)))
        if (present(parent_faces)) allocate(parent_faces(6*size(triangles, 2)))
        if (present(local_sectors)) allocate(local_sectors(6*size(triangles, 2)))
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
            refined_triangles(:, 6*face - 4) = [b, center, midpoint_ab]
            refined_triangles(:, 6*face - 3) = [b, midpoint_bc, center]
            refined_triangles(:, 6*face - 2) = [c, center, midpoint_bc]
            refined_triangles(:, 6*face - 1) = [c, midpoint_ca, center]
            refined_triangles(:, 6*face) = [a, center, midpoint_ca]
            if (present(parent_faces)) then
                parent_faces(6*face - 5:6*face) = face
            end if
            if (present(local_sectors)) then
                local_sectors(6*face - 5:6*face) = [1, 2, 3, 4, 5, 6]
            end if
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
