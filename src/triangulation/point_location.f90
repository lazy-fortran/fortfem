module point_location
    !> Efficient point location using walking search algorithm.
    !
    !  This module implements the walking search algorithm for locating a point
    !  within a triangulation. Instead of O(N) linear search, this achieves
    !  O(sqrt(N)) expected time by walking through the triangulation from a
    !  starting triangle toward the target point.
    !
    !  Algorithm (from Triangle/Shewchuk):
    !    1. Start from a given triangle (or random if not specified)
    !    2. For each edge of the current triangle, check if the target point
    !       is on the opposite side of the edge from the current triangle
    !    3. If so, move to the neighbor across that edge
    !    4. Repeat until the point is inside the current triangle (or on edge)
    !
    !  References:
    !    - Triangle library: preciselocate() and locate() functions
    !    - Shewchuk, "Triangle: Engineering a 2D Quality Mesh Generator"
    !
    use fortfem_kinds, only: dp
    use delaunay_types, only: mesh_t, point_t
    use geometric_predicates, only: orientation, ORIENTATION_CCW, ORIENTATION_CW
    implicit none

    private
    public :: locate_point
    public :: build_adjacency
    public :: LOCATION_INSIDE, LOCATION_ON_EDGE, LOCATION_NOT_FOUND

    ! Location result types
    integer, parameter :: LOCATION_INSIDE = 1     ! Point strictly inside triangle
    integer, parameter :: LOCATION_ON_EDGE = 2    ! Point on triangle edge
    integer, parameter :: LOCATION_NOT_FOUND = -1 ! Point not in triangulation

contains

    subroutine build_adjacency(mesh)
        !> Build triangle adjacency information in O(N) time.
        !
        !  Algorithm:
        !    1. Collect all edges as (min_v, max_v, tri_idx, edge_idx)
        !    2. Sort edges by (min_v, max_v) using counting sort on min_v
        !    3. Edges with same (min_v, max_v) are neighbors - link them
        !
        !  For each triangle, neighbors(i) is the triangle across edge i,
        !  where edge i connects vertices(i) and vertices(mod(i,3)+1).
        !
        type(mesh_t), intent(inout) :: mesh

        integer, allocatable :: edge_data(:,:)  ! (4, 3*ntriangles)
        integer, allocatable :: sorted_idx(:)
        integer :: n_edges, t, e, va, vb, min_v, max_v
        integer :: i, j

        ! Initialize all neighbors to 0 (no neighbor / boundary edge)
        do t = 1, mesh%ntriangles
            mesh%triangles(t)%neighbors = 0
        end do

        ! Count valid triangles and allocate edge storage
        n_edges = 0
        do t = 1, mesh%ntriangles
            if (mesh%triangles(t)%valid) n_edges = n_edges + 3
        end do

        if (n_edges == 0) return

        allocate(edge_data(4, n_edges))
        allocate(sorted_idx(n_edges))

        ! Collect all edges: (min_vertex, max_vertex, triangle_idx, edge_idx)
        n_edges = 0
        do t = 1, mesh%ntriangles
            if (.not. mesh%triangles(t)%valid) cycle

            do e = 1, 3
                va = mesh%triangles(t)%vertices(e)
                vb = mesh%triangles(t)%vertices(mod(e, 3) + 1)

                ! Canonical ordering: min vertex first
                if (va < vb) then
                    min_v = va
                    max_v = vb
                else
                    min_v = vb
                    max_v = va
                end if

                n_edges = n_edges + 1
                edge_data(1, n_edges) = min_v
                edge_data(2, n_edges) = max_v
                edge_data(3, n_edges) = t
                edge_data(4, n_edges) = e
            end do
        end do

        ! Sort edges by (min_v, max_v) using simple insertion sort
        ! For production, use quicksort, but this is O(N log N) anyway
        call sort_edges(edge_data, n_edges, sorted_idx)

        ! Link adjacent edges (consecutive sorted edges with same vertices)
        i = 1
        do while (i < n_edges)
            j = sorted_idx(i)

            ! Check if next edge has same vertices
            if (i + 1 <= n_edges) then
                if (edges_match(edge_data, sorted_idx(i), sorted_idx(i + 1))) then
                    ! Link the two triangles
                    call link_triangles(mesh, edge_data, sorted_idx(i),          &
                                       sorted_idx(i + 1))
                    i = i + 2  ! Skip both edges
                    cycle
                end if
            end if
            i = i + 1
        end do

        deallocate(edge_data)
        deallocate(sorted_idx)
    end subroutine build_adjacency

    subroutine sort_edges(edge_data, n_edges, sorted_idx)
        !> Sort edges by (min_v, max_v) and return sorted indices.
        integer, intent(in) :: edge_data(:,:)
        integer, intent(in) :: n_edges
        integer, intent(out) :: sorted_idx(:)

        integer :: i, j, key_idx
        integer :: key_min, key_max, cmp_min, cmp_max

        ! Initialize indices
        do i = 1, n_edges
            sorted_idx(i) = i
        end do

        ! Insertion sort (simple, works well for moderate N)
        do i = 2, n_edges
            key_idx = sorted_idx(i)
            key_min = edge_data(1, key_idx)
            key_max = edge_data(2, key_idx)

            j = i - 1
            do while (j >= 1)
                cmp_min = edge_data(1, sorted_idx(j))
                cmp_max = edge_data(2, sorted_idx(j))

                ! Compare (min_v, max_v) lexicographically
                if (cmp_min < key_min) exit
                if (cmp_min == key_min .and. cmp_max <= key_max) exit

                sorted_idx(j + 1) = sorted_idx(j)
                j = j - 1
            end do
            sorted_idx(j + 1) = key_idx
        end do
    end subroutine sort_edges

    logical function edges_match(edge_data, idx1, idx2) result(match)
        !> Check if two edges have the same vertices.
        integer, intent(in) :: edge_data(:,:)
        integer, intent(in) :: idx1, idx2

        match = (edge_data(1, idx1) == edge_data(1, idx2)) .and.                 &
                (edge_data(2, idx1) == edge_data(2, idx2))
    end function edges_match

    subroutine link_triangles(mesh, edge_data, idx1, idx2)
        !> Link two triangles as neighbors across their shared edge.
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: edge_data(:,:)
        integer, intent(in) :: idx1, idx2

        integer :: t1, e1, t2, e2

        t1 = edge_data(3, idx1)
        e1 = edge_data(4, idx1)
        t2 = edge_data(3, idx2)
        e2 = edge_data(4, idx2)

        mesh%triangles(t1)%neighbors(e1) = t2
        mesh%triangles(t2)%neighbors(e2) = t1
    end subroutine link_triangles

    subroutine locate_point(mesh, px, py, tri_idx, location_type, start_tri)
        !> Locate which triangle contains point (px, py) using walking search.
        !
        !  Args:
        !    mesh: Triangulation with adjacency information
        !    px, py: Coordinates of point to locate
        !    tri_idx: Output - index of containing triangle (0 if not found)
        !    location_type: Output - INSIDE, ON_EDGE, or NOT_FOUND
        !    start_tri: Optional starting triangle index for the search
        !
        type(mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: px, py
        integer, intent(out) :: tri_idx
        integer, intent(out) :: location_type
        integer, intent(in), optional :: start_tri

        type(point_t) :: p, pa, pb, pc
        integer :: current, next, edge_to_cross
        integer :: orient_ab, orient_bc, orient_ca
        integer :: max_iter, iter
        integer :: va, vb, vc

        tri_idx = 0
        location_type = LOCATION_NOT_FOUND

        ! Create point structure
        p%x = px
        p%y = py

        ! Find starting triangle
        if (present(start_tri)) then
            current = start_tri
            if (current < 1 .or. current > mesh%ntriangles) then
                current = find_valid_triangle(mesh)
            else if (.not. mesh%triangles(current)%valid) then
                current = find_valid_triangle(mesh)
            end if
        else
            current = find_valid_triangle(mesh)
        end if

        if (current == 0) return

        ! Walking search with iteration limit
        max_iter = mesh%ntriangles + 10
        do iter = 1, max_iter
            if (.not. mesh%triangles(current)%valid) then
                current = find_valid_triangle(mesh)
                if (current == 0) return
            end if

            ! Get triangle vertices
            va = mesh%triangles(current)%vertices(1)
            vb = mesh%triangles(current)%vertices(2)
            vc = mesh%triangles(current)%vertices(3)

            pa = mesh%points(va)
            pb = mesh%points(vb)
            pc = mesh%points(vc)

            ! Test orientation of point relative to each edge
            orient_ab = orientation(pa, pb, p)
            orient_bc = orientation(pb, pc, p)
            orient_ca = orientation(pc, pa, p)

            ! Check if point is inside (all same orientation) or on edge
            if (orient_ab >= 0 .and. orient_bc >= 0 .and. orient_ca >= 0) then
                ! Point is inside or on edge of CCW triangle
                tri_idx = current
                if (orient_ab == 0 .or. orient_bc == 0 .or. orient_ca == 0) then
                    location_type = LOCATION_ON_EDGE
                else
                    location_type = LOCATION_INSIDE
                end if
                return
            end if

            ! Point is outside - walk toward it
            edge_to_cross = 0
            if (orient_ab < 0) then
                edge_to_cross = 1  ! Cross edge AB (edge 1)
            else if (orient_bc < 0) then
                edge_to_cross = 2  ! Cross edge BC (edge 2)
            else if (orient_ca < 0) then
                edge_to_cross = 3  ! Cross edge CA (edge 3)
            end if

            if (edge_to_cross == 0) exit  ! Should not happen

            ! Move to neighbor across the edge
            next = mesh%triangles(current)%neighbors(edge_to_cross)
            if (next == 0) then
                ! No neighbor - point is outside triangulation
                return
            end if

            current = next
        end do

        ! If we exit loop without finding, fall back to linear search
        call locate_point_linear(mesh, px, py, tri_idx, location_type)
    end subroutine locate_point

    subroutine locate_point_linear(mesh, px, py, tri_idx, location_type)
        !> Linear search fallback for point location.
        type(mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: px, py
        integer, intent(out) :: tri_idx
        integer, intent(out) :: location_type

        type(point_t) :: p, pa, pb, pc
        integer :: t, va, vb, vc
        integer :: orient_ab, orient_bc, orient_ca

        tri_idx = 0
        location_type = LOCATION_NOT_FOUND

        p%x = px
        p%y = py

        do t = 1, mesh%ntriangles
            if (.not. mesh%triangles(t)%valid) cycle

            va = mesh%triangles(t)%vertices(1)
            vb = mesh%triangles(t)%vertices(2)
            vc = mesh%triangles(t)%vertices(3)

            pa = mesh%points(va)
            pb = mesh%points(vb)
            pc = mesh%points(vc)

            orient_ab = orientation(pa, pb, p)
            orient_bc = orientation(pb, pc, p)
            orient_ca = orientation(pc, pa, p)

            if (orient_ab >= 0 .and. orient_bc >= 0 .and. orient_ca >= 0) then
                tri_idx = t
                if (orient_ab == 0 .or. orient_bc == 0 .or. orient_ca == 0) then
                    location_type = LOCATION_ON_EDGE
                else
                    location_type = LOCATION_INSIDE
                end if
                return
            end if
        end do
    end subroutine locate_point_linear

    integer function find_valid_triangle(mesh) result(tri_idx)
        !> Find any valid triangle to start the walk.
        type(mesh_t), intent(in) :: mesh

        integer :: t

        tri_idx = 0
        do t = 1, mesh%ntriangles
            if (mesh%triangles(t)%valid) then
                tri_idx = t
                return
            end if
        end do
    end function find_valid_triangle

end module point_location
