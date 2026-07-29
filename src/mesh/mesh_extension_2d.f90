module fortfem_mesh_extension_2d
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use triangulation_fortran, only: cleanup_triangulation, &
        triangulate_with_hole_fortran, triangulation_result_t
    implicit none
    private

    public :: extend_triangle_mesh

contains

    subroutine extend_triangle_mesh(core, outer_vertices, extended, status)
        type(mesh_2d_t), intent(inout) :: core
        real(dp), intent(in) :: outer_vertices(:, :)
        type(mesh_2d_t), intent(out) :: extended
        integer, intent(out) :: status

        type(triangulation_result_t) :: annulus
        integer, allocatable :: boundary_vertices(:), point_map(:)
        integer, allocatable :: segments(:, :)
        real(dp), allocatable :: annulus_points(:, :), hole_point(:, :)
        integer :: annulus_status, inner_count, new_count, outer_count
        integer :: point, triangle

        status = -1
        if (.not. valid_core_mesh(core)) return
        if (size(outer_vertices, 1) /= 2) return
        outer_count = size(outer_vertices, 2)
        if (outer_count < 3) return

        if (.not. allocated(core%edges)) call core%build_edge_connectivity()
        if (.not. allocated(core%boundary_edges)) call core%find_boundary()
        inner_count = core%n_boundary_edges
        if (inner_count < 3) return
        call order_boundary_vertices(core, boundary_vertices, status)
        if (status /= 0) return

        allocate(annulus_points(2, inner_count + outer_count))
        do point = 1, inner_count
            annulus_points(:, point) = &
                core%vertices(:, boundary_vertices(point))
        end do
        annulus_points(:, inner_count + 1:) = outer_vertices
        call make_loop_segments(inner_count, outer_count, segments)
        allocate(hole_point(2, 1))
        hole_point(:, 1) = sum( &
            core%vertices(:, core%triangles(:, 1)), dim=2) / 3.0_dp

        call triangulate_with_hole_fortran(annulus_points, segments, &
            hole_point, annulus, annulus_status)
        if (annulus_status /= 0 .or. annulus%ntriangles < 1) then
            status = -2
            call cleanup_triangulation(annulus)
            return
        end if

        allocate(point_map(annulus%npoints))
        call map_annulus_points( &
            core, boundary_vertices, annulus, point_map, new_count, status)
        if (status /= 0) then
            call cleanup_triangulation(annulus)
            return
        end if

        extended%n_vertices = core%n_vertices + new_count
        extended%n_triangles = core%n_triangles + annulus%ntriangles
        extended%has_triangles = .true.
        allocate(extended%vertices(2, extended%n_vertices))
        allocate(extended%triangles(3, extended%n_triangles))
        extended%vertices(:, 1:core%n_vertices) = core%vertices
        extended%triangles(:, 1:core%n_triangles) = core%triangles
        do point = 1, annulus%npoints
            if (point_map(point) > core%n_vertices) then
                extended%vertices(:, point_map(point)) = annulus%points(:, point)
            end if
        end do
        do triangle = 1, annulus%ntriangles
            extended%triangles(:, core%n_triangles + triangle) = &
                point_map(annulus%triangles(:, triangle))
            call orient_triangle_counterclockwise( &
                extended, core%n_triangles + triangle)
        end do
        call extended%build_edge_connectivity()
        call extended%build_edge_dof_numbering()
        call cleanup_triangulation(annulus)
        status = 0
    end subroutine extend_triangle_mesh

    logical function valid_core_mesh(core) result(valid)
        type(mesh_2d_t), intent(in) :: core

        real(dp) :: determinant
        integer :: triangle

        valid = core%has_triangles .and. core%n_vertices >= 3 .and. &
            core%n_triangles >= 1
        if (.not. valid) return
        if (.not. allocated(core%vertices) .or. &
            .not. allocated(core%triangles)) then
            valid = .false.
            return
        end if
        do triangle = 1, core%n_triangles
            determinant = triangle_determinant(core, triangle)
            if (determinant <= 0.0_dp) then
                valid = .false.
                return
            end if
        end do
    end function valid_core_mesh

    subroutine order_boundary_vertices(mesh, vertices, status)
        type(mesh_2d_t), intent(in) :: mesh
        integer, allocatable, intent(out) :: vertices(:)
        integer, intent(out) :: status

        logical, allocatable :: used(:)
        integer :: boundary_edge, edge, position, current, next_vertex

        allocate(vertices(mesh%n_boundary_edges))
        allocate(used(mesh%n_boundary_edges))
        used = .false.
        boundary_edge = mesh%boundary_edges(1)
        vertices(1) = mesh%edges(1, boundary_edge)
        vertices(2) = mesh%edges(2, boundary_edge)
        used(1) = .true.

        do position = 3, mesh%n_boundary_edges
            current = vertices(position - 1)
            next_vertex = 0
            do edge = 1, mesh%n_boundary_edges
                if (used(edge)) cycle
                boundary_edge = mesh%boundary_edges(edge)
                if (mesh%edges(1, boundary_edge) == current) then
                    next_vertex = mesh%edges(2, boundary_edge)
                else if (mesh%edges(2, boundary_edge) == current) then
                    next_vertex = mesh%edges(1, boundary_edge)
                else
                    cycle
                end if
                used(edge) = .true.
                exit
            end do
            if (next_vertex == 0) then
                status = -1
                return
            end if
            vertices(position) = next_vertex
        end do
        status = -1
        current = vertices(mesh%n_boundary_edges)
        do edge = 1, mesh%n_boundary_edges
            if (used(edge)) cycle
            boundary_edge = mesh%boundary_edges(edge)
            if (all(mesh%edges(:, boundary_edge) == &
                [min(current, vertices(1)), max(current, vertices(1))])) then
                status = 0
                return
            end if
        end do
    end subroutine order_boundary_vertices

    subroutine make_loop_segments(inner_count, outer_count, segments)
        integer, intent(in) :: inner_count, outer_count
        integer, allocatable, intent(out) :: segments(:, :)

        integer :: point

        allocate(segments(2, inner_count + outer_count))
        do point = 1, inner_count
            segments(:, point) = [point, mod(point, inner_count) + 1]
        end do
        do point = 1, outer_count
            segments(:, inner_count + point) = [ &
                inner_count + point, &
                inner_count + mod(point, outer_count) + 1]
        end do
    end subroutine make_loop_segments

    subroutine map_annulus_points( &
            core, boundary_vertices, annulus, point_map, new_count, status)
        type(mesh_2d_t), intent(in) :: core
        integer, intent(in) :: boundary_vertices(:)
        type(triangulation_result_t), intent(in) :: annulus
        integer, intent(out) :: point_map(:), new_count, status

        integer :: boundary_point, point
        logical :: matched

        new_count = 0
        do point = 1, annulus%npoints
            matched = .false.
            do boundary_point = 1, size(boundary_vertices)
                if (same_point(annulus%points(:, point), &
                    core%vertices(:, boundary_vertices(boundary_point)))) then
                    point_map(point) = boundary_vertices(boundary_point)
                    matched = .true.
                    exit
                end if
            end do
            if (.not. matched) then
                new_count = new_count + 1
                point_map(point) = core%n_vertices + new_count
            end if
        end do
        status = 0
        if (count(point_map <= core%n_vertices) /= size(boundary_vertices)) then
            status = -3
        end if
    end subroutine map_annulus_points

    pure logical function same_point(a, b) result(same)
        real(dp), intent(in) :: a(2), b(2)

        real(dp) :: scale

        scale = max(1.0_dp, maxval(abs(a)), maxval(abs(b)))
        same = maxval(abs(a - b)) <= 32.0_dp * epsilon(1.0_dp) * scale
    end function same_point

    subroutine orient_triangle_counterclockwise(mesh, triangle)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: triangle

        integer :: temporary

        if (triangle_determinant(mesh, triangle) < 0.0_dp) then
            temporary = mesh%triangles(2, triangle)
            mesh%triangles(2, triangle) = mesh%triangles(3, triangle)
            mesh%triangles(3, triangle) = temporary
        end if
    end subroutine orient_triangle_counterclockwise

    pure real(dp) function triangle_determinant(mesh, triangle) result(value)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle

        value = &
            (mesh%vertices(1, mesh%triangles(2, triangle)) - &
            mesh%vertices(1, mesh%triangles(1, triangle))) * &
            (mesh%vertices(2, mesh%triangles(3, triangle)) - &
            mesh%vertices(2, mesh%triangles(1, triangle))) - &
            (mesh%vertices(1, mesh%triangles(3, triangle)) - &
            mesh%vertices(1, mesh%triangles(1, triangle))) * &
            (mesh%vertices(2, mesh%triangles(2, triangle)) - &
            mesh%vertices(2, mesh%triangles(1, triangle)))
    end function triangle_determinant

end module fortfem_mesh_extension_2d
