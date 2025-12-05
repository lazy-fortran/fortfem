module bowyer_watson
    !> Bowyer-Watson algorithm for Delaunay triangulation.
    !
    !  This implementation uses:
    !    - Robust geometric predicates (integer coordinate strategy)
    !    - O(sqrt N) point location via walking search
    !    - BFS cavity expansion from seed triangle (O(1) average cavity size)
    !    - Incremental adjacency maintenance
    !
    !  Overall complexity: O(N log N) expected for N point insertions.
    !
    !  References:
    !    - Bowyer, A. (1981). Computing Dirichlet tessellations.
    !    - Watson, D. F. (1981). Computing the n-dimensional Delaunay tessellation.
    !
    use fortfem_kinds, only: dp
    use delaunay_types
    use geometric_predicates
    use point_location, only: build_adjacency, locate_point, LOCATION_INSIDE,    &
        LOCATION_ON_EDGE, LOCATION_NOT_FOUND
    implicit none

    private
    public :: delaunay_triangulate, insert_point
    public :: create_super_triangle, remove_super_triangle
    public :: find_cavity, fill_cavity

    ! Maximum expected cavity size (for static allocation)
    integer, parameter :: MAX_CAVITY_SIZE = 64

    ! Last triangle used for point location (for walking search warm start)
    integer, save :: last_triangle = 1

contains

    subroutine delaunay_triangulate(input_points, mesh)
        !> Main Delaunay triangulation routine using Bowyer-Watson algorithm.
        !  Robust predicates are always used.
        real(dp), intent(in) :: input_points(:,:)  ! (2, npoints)
        type(mesh_t), intent(out) :: mesh

        integer :: i, npoints
        integer :: super_idx, point_idx

        npoints = size(input_points, 2)

        ! Initialize robust predicates for this point set
        call init_robust_predicates(input_points, npoints)

        ! Initialize mesh with appropriate size
        call create_mesh(mesh, npoints + 3, 6 * npoints, 3 * npoints)

        ! Create super-triangle containing all points
        call create_super_triangle(input_points, mesh, super_idx)

        ! Reset walking search state
        last_triangle = 1

        ! Add input points to mesh
        do i = 1, npoints
            point_idx = add_point(mesh, input_points(1, i), input_points(2, i), i)
        end do

        ! Build initial adjacency for super-triangle once
        call build_adjacency(mesh)

        ! Insert each point using Bowyer-Watson algorithm
        ! Adjacency is maintained incrementally in fill_cavity (O(1) overhead)
        do i = 4, mesh%npoints  ! Start after super-triangle vertices
            ! write(*,*) "Inserting point ", i
            call insert_point(mesh, i)
        end do

        ! Remove super-triangle and its associated triangles
        call remove_super_triangle(mesh)
    end subroutine delaunay_triangulate

    subroutine create_super_triangle(input_points, mesh, super_tri_idx)
        !> Create a super-triangle that contains all input points.
        real(dp), intent(in) :: input_points(:,:)
        type(mesh_t), intent(inout) :: mesh
        integer, intent(out) :: super_tri_idx

        real(dp) :: min_x, max_x, min_y, max_y
        real(dp) :: dx, dy, delta_max, x_mid, y_mid
        integer :: p1, p2, p3

        ! Find bounding box of all points
        min_x = minval(input_points(1, :))
        max_x = maxval(input_points(1, :))
        min_y = minval(input_points(2, :))
        max_y = maxval(input_points(2, :))

        dx = max_x - min_x
        dy = max_y - min_y
        delta_max = max(dx, dy)
        if (delta_max < 1.0e-14_dp) delta_max = 1.0_dp

        x_mid = (min_x + max_x) / 2.0_dp
        y_mid = (min_y + max_y) / 2.0_dp

        ! Create super-triangle vertices (much larger than bounding box)
        p1 = add_point(mesh, x_mid - 20.0_dp * delta_max, y_mid - delta_max, -1)
        p2 = add_point(mesh, x_mid + 20.0_dp * delta_max, y_mid - delta_max, -2)
        p3 = add_point(mesh, x_mid, y_mid + 20.0_dp * delta_max, -3)

        ! Store super-triangle vertex indices
        mesh%super_vertices(1) = p1
        mesh%super_vertices(2) = p2
        mesh%super_vertices(3) = p3

        ! Create super-triangle
        super_tri_idx = add_triangle(mesh, p1, p2, p3)
    end subroutine create_super_triangle

    subroutine remove_super_triangle(mesh)
        !> Remove super-triangle vertices and all triangles containing them.
        type(mesh_t), intent(inout) :: mesh

        integer :: i, j, v, valid_triangle_count
        logical :: contains_super_vertex

        valid_triangle_count = 0

        ! Mark triangles containing super-triangle vertices as invalid
        do i = 1, mesh%ntriangles
            if (.not. mesh%triangles(i)%valid) cycle

            contains_super_vertex = .false.
            do j = 1, 3
                v = mesh%triangles(i)%vertices(j)
                if (v == mesh%super_vertices(1) .or.                             &
                    v == mesh%super_vertices(2) .or.                             &
                    v == mesh%super_vertices(3)) then
                    contains_super_vertex = .true.
                    exit
                end if
            end do

            if (contains_super_vertex) then
                mesh%triangles(i)%valid = .false.
            else
                valid_triangle_count = valid_triangle_count + 1
            end if
        end do

        ! If no valid triangles remain, re-triangulate real vertices.
        ! This should only happen for degenerate cases (all points collinear, etc.)
        if (valid_triangle_count == 0) then
            call create_triangles_from_real_vertices(mesh)
        end if

        ! Mark super-triangle vertices as invalid
        do i = 1, 3
            if (mesh%super_vertices(i) > 0 .and.                                 &
                mesh%super_vertices(i) <= mesh%npoints) then
                mesh%points(mesh%super_vertices(i))%valid = .false.
            end if
        end do
    end subroutine remove_super_triangle

    subroutine create_triangles_from_real_vertices(mesh)
        !> Create triangles using only real vertices (handles N >= 3).
        !
        !  For convex point sets, this creates a fan triangulation from
        !  the convex hull. For N=3, creates a single triangle.
        !
        type(mesh_t), intent(inout) :: mesh

        integer, allocatable :: real_verts(:)
        integer, allocatable :: hull(:)
        integer :: n_real, n_hull, i, tri_idx
        type(point_t) :: p1, p2, p3
        integer :: orient

        ! Collect all real vertices (positive IDs)
        allocate(real_verts(mesh%npoints))
        n_real = 0

        do i = 1, mesh%npoints
            if (mesh%points(i)%valid .and. mesh%points(i)%id > 0) then
                n_real = n_real + 1
                real_verts(n_real) = i
            end if
        end do

        if (n_real < 3) then
            deallocate(real_verts)
            return
        end if

        ! Compute convex hull
        call convex_hull(mesh, real_verts, n_real, hull, n_hull)

        if (n_hull < 3) then
            deallocate(real_verts)
            if (allocated(hull)) deallocate(hull)
            return
        end if

        ! Fan triangulation from first hull vertex
        do i = 2, n_hull - 1
            p1 = mesh%points(hull(1))
            p2 = mesh%points(hull(i))
            p3 = mesh%points(hull(i + 1))

            orient = orientation(p1, p2, p3)
            if (orient == ORIENTATION_CCW) then
                tri_idx = add_triangle(mesh, hull(1), hull(i), hull(i + 1))
            else if (orient == ORIENTATION_CW) then
                tri_idx = add_triangle(mesh, hull(1), hull(i + 1), hull(i))
            end if
        end do

        deallocate(real_verts)
        deallocate(hull)
    end subroutine create_triangles_from_real_vertices

    subroutine convex_hull(mesh, verts, n_verts, hull, n_hull)
        !> Compute convex hull using gift wrapping algorithm.
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: verts(:)
        integer, intent(in) :: n_verts
        integer, allocatable, intent(out) :: hull(:)
        integer, intent(out) :: n_hull

        integer :: i, start_idx, current, next, candidate
        integer :: orient
        type(point_t) :: pc, pn, pcand
        logical, allocatable :: on_hull(:)

        allocate(hull(n_verts))
        allocate(on_hull(n_verts))
        on_hull = .false.
        n_hull = 0

        ! Find leftmost point
        start_idx = 1
        do i = 2, n_verts
            if (mesh%points(verts(i))%x < mesh%points(verts(start_idx))%x) then
                start_idx = i
            else if (mesh%points(verts(i))%x == mesh%points(verts(start_idx))%x   &
                .and. mesh%points(verts(i))%y < mesh%points(verts(start_idx))%y) then
                start_idx = i
            end if
        end do

        current = start_idx
        do
            n_hull = n_hull + 1
            hull(n_hull) = verts(current)
            on_hull(current) = .true.

            ! Find next hull vertex (most counter-clockwise from current)
            next = 1
            if (next == current) next = 2
            if (next > n_verts) exit

            pc = mesh%points(verts(current))
            pn = mesh%points(verts(next))

            do candidate = 1, n_verts
                if (candidate == current .or. candidate == next) cycle

                pcand = mesh%points(verts(candidate))
                orient = orientation(pc, pn, pcand)

                if (orient == ORIENTATION_CCW) then
                    next = candidate
                    pn = mesh%points(verts(next))
                end if
            end do

            current = next
            if (current == start_idx) exit
            if (n_hull >= n_verts) exit
        end do

        deallocate(on_hull)
    end subroutine convex_hull

    subroutine insert_point(mesh, point_idx)
        !> Insert a point into existing triangulation using Bowyer-Watson.
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: point_idx

        integer :: cavity_triangles(MAX_CAVITY_SIZE)
        integer :: cavity_edges(2, MAX_CAVITY_SIZE * 3)
        integer :: external_neighbors(MAX_CAVITY_SIZE * 3)
        integer :: ncavity_triangles, ncavity_edges

        ! Find triangles whose circumcircles contain the new point
        call find_cavity(mesh, point_idx, cavity_triangles, ncavity_triangles)

        if (ncavity_triangles == -1) return  ! Handled via edge splitting

        if (ncavity_triangles == 0) return

        ! Find the boundary of the cavity and external neighbors
        call find_cavity_boundary(mesh, cavity_triangles, ncavity_triangles,     &
                                  cavity_edges, external_neighbors, ncavity_edges)

        ! Remove triangles in the cavity
        call remove_cavity_triangles(mesh, cavity_triangles, ncavity_triangles)

        ! Create new triangles with incremental adjacency updates
        call fill_cavity(mesh, point_idx, cavity_edges, external_neighbors,      &
                        ncavity_edges)
    end subroutine insert_point

    subroutine find_cavity(mesh, point_idx, cavity_triangles, ncavity_triangles)
        !> Find all triangles whose circumcircles contain the given point.
        !
        !  Uses walking search to find seed triangle, then BFS expansion.
        !  This achieves O(sqrt N) for point location + O(k) for cavity
        !  where k is cavity size (typically constant).
        !
        !  When the point lies on an edge (LOCATION_ON_EDGE), we directly split
        !  the edge instead of using cavity-based insertion. This avoids creating
        !  degenerate triangles with collinear points.
        !
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: point_idx
        integer, intent(out) :: cavity_triangles(MAX_CAVITY_SIZE)
        integer, intent(out) :: ncavity_triangles

        type(point_t) :: point, pa, pb, pc
        integer :: seed_tri, loc_type
        integer :: queue(MAX_CAVITY_SIZE)
        logical :: visited(mesh%ntriangles)
        integer :: queue_head, queue_tail
        integer :: t, neighbor, e
        integer :: edge_neighbor

        point = mesh%points(point_idx)
        ncavity_triangles = 0
        visited = .false.

        ! Step 1: Find seed triangle using walking search
        ! locate_point automatically handles linear fallback if walk fails.
        call locate_point(mesh, point%x, point%y, seed_tri, loc_type,            &
                         last_triangle)

        if (seed_tri == 0) return

        ! Update last_triangle for next search warm start
        last_triangle = seed_tri

        ! Special handling for points on edges (segment midpoints)
        ! Use direct edge splitting instead of cavity-based insertion to avoid
        ! creating degenerate triangles with collinear points.
        if (loc_type == LOCATION_ON_EDGE) then
            call split_edge_triangles(mesh, seed_tri, point_idx)
            ! Return -1 to indicate handled via edge splitting
            ncavity_triangles = -1
            return
        end if

        ! Step 2: BFS expansion from seed triangle (standard Bowyer-Watson)
        queue_head = 1
        queue_tail = 1
        queue(queue_tail) = seed_tri
        visited(seed_tri) = .true.

        do while (queue_head <= queue_tail)
            t = queue(queue_head)
            queue_head = queue_head + 1

            if (.not. mesh%triangles(t)%valid) cycle

            ! Check if point is in circumcircle of triangle t
            pa = mesh%points(mesh%triangles(t)%vertices(1))
            pb = mesh%points(mesh%triangles(t)%vertices(2))
            pc = mesh%points(mesh%triangles(t)%vertices(3))

            if (in_circle(pa, pb, pc, point)) then
                ! Add to cavity
                if (ncavity_triangles < MAX_CAVITY_SIZE) then
                    ncavity_triangles = ncavity_triangles + 1
                    cavity_triangles(ncavity_triangles) = t
                end if

                ! Expand to neighbors
                do e = 1, 3
                    neighbor = mesh%triangles(t)%neighbors(e)
                    if (neighbor > 0 .and. neighbor <= mesh%ntriangles) then
                        if (.not. visited(neighbor)) then
                            visited(neighbor) = .true.
                            if (queue_tail < MAX_CAVITY_SIZE) then
                                queue_tail = queue_tail + 1
                                queue(queue_tail) = neighbor
                            end if
                        end if
                    end if
                end do
            end if
        end do

        ! Failsafe: If point is strictly inside a triangle, that triangle MUST be
        ! in the cavity (geometry guarantee). If robust predicates failed to
        ! identify it (e.g. point very close to edge), force it in.
        if (ncavity_triangles == 0 .and. loc_type == LOCATION_INSIDE) then
            ncavity_triangles = 1
            cavity_triangles(1) = seed_tri
        end if
    end subroutine find_cavity

    subroutine split_edge_triangles(mesh, seed_tri, point_idx)
        !> Split the edge containing a point by directly modifying triangles.
        !
        !  When a point lies exactly on an edge, we can't use the standard
        !  Bowyer-Watson cavity approach because it may create degenerate
        !  triangles. Instead, we directly split the edge:
        !
        !  Before:
        !       v3                v3
        !      /  \              /||\
        !     / t1 \            / || \
        !    /      \   =>     /  ||  \
        !   v1------v2        v1--p--v2
        !    \      /          \  ||  /
        !     \ t2 /            \ || /
        !      \  /              \||/
        !       v4                v4
        !
        !  After splitting, we have 4 triangles instead of 2.
        !
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: seed_tri
        integer, intent(in) :: point_idx

        type(point_t) :: point, p1, p2
        integer :: e, v1, v2, v3, v4, neighbor
        integer :: t1, t2, t3, t4
        integer :: edge_idx, neighbor_edge_idx
        real(dp) :: dx, dy, len2, px, py
        logical :: found_edge, has_neighbor

        point = mesh%points(point_idx)
        found_edge = .false.

        ! Find which edge the point lies on
        do e = 1, 3
            v1 = mesh%triangles(seed_tri)%vertices(e)
            v2 = mesh%triangles(seed_tri)%vertices(mod(e, 3) + 1)
            v3 = mesh%triangles(seed_tri)%vertices(mod(e + 1, 3) + 1)

            p1 = mesh%points(v1)
            p2 = mesh%points(v2)

            dx = p2%x - p1%x
            dy = p2%y - p1%y
            len2 = dx*dx + dy*dy

            if (len2 < 1.0e-20_dp) cycle

            px = point%x - p1%x
            py = point%y - p1%y

            ! Check perpendicular distance to edge
            if (abs(px*dy - py*dx) / sqrt(len2) < 1.0e-10_dp) then
                ! Point is on this edge
                edge_idx = e
                found_edge = .true.
                exit
            end if
        end do

        if (.not. found_edge) return

        ! Get the neighbor triangle across this edge
        neighbor = mesh%triangles(seed_tri)%neighbors(edge_idx)
        has_neighbor = (neighbor > 0 .and. neighbor <= mesh%ntriangles)
        if (has_neighbor) has_neighbor = mesh%triangles(neighbor)%valid

        ! Find v4 (the opposite vertex in the neighbor triangle)
        v4 = 0
        if (has_neighbor) then
            do e = 1, 3
                if (mesh%triangles(neighbor)%vertices(e) /= v1 .and.             &
                    mesh%triangles(neighbor)%vertices(e) /= v2) then
                    v4 = mesh%triangles(neighbor)%vertices(e)
                    neighbor_edge_idx = e
                    exit
                end if
            end do
        end if

        ! Mark original triangles as invalid
        mesh%triangles(seed_tri)%valid = .false.
        if (has_neighbor) mesh%triangles(neighbor)%valid = .false.

        ! Create new triangles
        ! Triangle 1: (v1, point, v3) - part of original seed_tri
        t1 = add_triangle(mesh, v1, point_idx, v3)

        ! Triangle 2: (point, v2, v3) - part of original seed_tri
        t2 = add_triangle(mesh, point_idx, v2, v3)

        if (has_neighbor .and. v4 > 0) then
            ! Triangle 3: (v2, point, v4) - part of original neighbor
            t3 = add_triangle(mesh, v2, point_idx, v4)

            ! Triangle 4: (point, v1, v4) - part of original neighbor
            t4 = add_triangle(mesh, point_idx, v1, v4)
        else
            t3 = 0
            t4 = 0
        end if

        ! Set up adjacencies for the new triangles
        ! t1: (v1, point, v3)
        !   Edge 1 (v1-point): adjacent to t4 if exists, else boundary
        !   Edge 2 (point-v3): adjacent to t2
        !   Edge 3 (v3-v1): adjacent to old neighbor of seed_tri edge (v3-v1)
        mesh%triangles(t1)%neighbors(1) = t4
        mesh%triangles(t1)%neighbors(2) = t2
        mesh%triangles(t1)%neighbors(3) = find_old_neighbor(mesh, seed_tri,      &
                                                            v3, v1, edge_idx)

        ! t2: (point, v2, v3)
        !   Edge 1 (point-v2): adjacent to t3 if exists, else boundary
        !   Edge 2 (v2-v3): adjacent to old neighbor of seed_tri edge (v2-v3)
        !   Edge 3 (v3-point): adjacent to t1
        mesh%triangles(t2)%neighbors(1) = t3
        mesh%triangles(t2)%neighbors(2) = find_old_neighbor(mesh, seed_tri,      &
                                                            v2, v3, edge_idx)
        mesh%triangles(t2)%neighbors(3) = t1

        if (has_neighbor .and. v4 > 0) then
            ! t3: (v2, point, v4)
            !   Edge 1 (v2-point): adjacent to t2
            !   Edge 2 (point-v4): adjacent to t4
            !   Edge 3 (v4-v2): adjacent to old neighbor of neighbor edge
            mesh%triangles(t3)%neighbors(1) = t2
            mesh%triangles(t3)%neighbors(2) = t4
            mesh%triangles(t3)%neighbors(3) = find_old_neighbor(mesh, neighbor,  &
                                                    v4, v2, neighbor_edge_idx)

            ! t4: (point, v1, v4)
            !   Edge 1 (point-v1): adjacent to t1
            !   Edge 2 (v1-v4): adjacent to old neighbor of neighbor edge
            !   Edge 3 (v4-point): adjacent to t3
            mesh%triangles(t4)%neighbors(1) = t1
            mesh%triangles(t4)%neighbors(2) = find_old_neighbor(mesh, neighbor,  &
                                                    v1, v4, neighbor_edge_idx)
            mesh%triangles(t4)%neighbors(3) = t3
        end if

        ! Update external neighbors to point to new triangles
        call update_external_neighbors_after_split(mesh, seed_tri, t1, t2,       &
                                                   v1, v2, v3, edge_idx)
        if (has_neighbor .and. v4 > 0) then
            call update_external_neighbors_after_split(mesh, neighbor, t3, t4,   &
                                                       v2, v1, v4,               &
                                                       neighbor_edge_idx)
        end if

        ! Update last_triangle for warm start
        last_triangle = t1
    end subroutine split_edge_triangles

    integer function find_old_neighbor(mesh, tri_idx, va, vb, skip_edge)
        !> Find the neighbor of tri_idx across edge (va, vb), skipping skip_edge.
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: tri_idx, va, vb, skip_edge

        integer :: e, ev1, ev2

        find_old_neighbor = 0

        do e = 1, 3
            if (e == skip_edge) cycle
            ev1 = mesh%triangles(tri_idx)%vertices(e)
            ev2 = mesh%triangles(tri_idx)%vertices(mod(e, 3) + 1)

            if ((ev1 == va .and. ev2 == vb) .or.                                 &
                (ev1 == vb .and. ev2 == va)) then
                find_old_neighbor = mesh%triangles(tri_idx)%neighbors(e)
                return
            end if
        end do
    end function find_old_neighbor

    subroutine update_external_neighbors_after_split(mesh, old_tri, new_t1,      &
                                                     new_t2, v1, v2, v3,         &
                                                     split_edge)
        !> Update external triangles to point to new triangles after edge split.
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: old_tri, new_t1, new_t2, v1, v2, v3, split_edge

        integer :: e, neighbor, ev1, ev2

        do e = 1, 3
            if (e == split_edge) cycle

            neighbor = mesh%triangles(old_tri)%neighbors(e)
            if (neighbor <= 0 .or. neighbor > mesh%ntriangles) cycle
            if (.not. mesh%triangles(neighbor)%valid) cycle

            ev1 = mesh%triangles(old_tri)%vertices(e)
            ev2 = mesh%triangles(old_tri)%vertices(mod(e, 3) + 1)

            ! Determine which new triangle inherits this edge
            if ((ev1 == v3 .and. ev2 == v1) .or. (ev1 == v1 .and. ev2 == v3)) then
                call update_neighbor_link(mesh, neighbor, ev1, ev2, new_t1)
            else if ((ev1 == v2 .and. ev2 == v3) .or.                            &
                     (ev1 == v3 .and. ev2 == v2)) then
                call update_neighbor_link(mesh, neighbor, ev1, ev2, new_t2)
            end if
        end do
    end subroutine update_external_neighbors_after_split

    subroutine find_cavity_boundary(mesh, cavity_triangles, ncavity_triangles,   &
                                    cavity_edges, external_neighbors,            &
                                    ncavity_edges)
        !> Find boundary edges of the cavity and their external neighbors.
        !
        !  For each boundary edge, we also store which triangle is on the
        !  outside of the cavity (needed for incremental adjacency updates).
        !
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: cavity_triangles(MAX_CAVITY_SIZE)
        integer, intent(in) :: ncavity_triangles
        integer, intent(out) :: cavity_edges(2, MAX_CAVITY_SIZE * 3)
        integer, intent(out) :: external_neighbors(MAX_CAVITY_SIZE * 3)
        integer, intent(out) :: ncavity_edges

        integer :: all_edges(2, MAX_CAVITY_SIZE * 3)
        integer :: edge_tri(MAX_CAVITY_SIZE * 3)      ! Which triangle owns edge
        integer :: edge_idx_in_tri(MAX_CAVITY_SIZE * 3)  ! Edge index in triangle
        integer :: edge_count(MAX_CAVITY_SIZE * 3)
        integer :: nedges_total, i, j, t, v1, v2
        integer :: edge_idx, tmp, neighbor
        logical :: found, in_cavity

        nedges_total = 0

        do i = 1, ncavity_triangles
            t = cavity_triangles(i)
            do j = 1, 3
                v1 = mesh%triangles(t)%vertices(j)
                v2 = mesh%triangles(t)%vertices(mod(j, 3) + 1)

                ! Canonical ordering
                if (v1 > v2) then
                    tmp = v1
                    v1 = v2
                    v2 = tmp
                end if

                ! Check if edge already exists
                found = .false.
                do edge_idx = 1, nedges_total
                    if (all_edges(1, edge_idx) == v1 .and.                       &
                        all_edges(2, edge_idx) == v2) then
                        edge_count(edge_idx) = edge_count(edge_idx) + 1
                        found = .true.
                        exit
                    end if
                end do

                if (.not. found) then
                    nedges_total = nedges_total + 1
                    all_edges(1, nedges_total) = v1
                    all_edges(2, nedges_total) = v2
                    edge_tri(nedges_total) = t
                    edge_idx_in_tri(nedges_total) = j
                    edge_count(nedges_total) = 1
                end if
            end do
        end do

        ! Collect boundary edges and their external neighbors
        ncavity_edges = 0
        do i = 1, nedges_total
            if (edge_count(i) == 1) then
                ncavity_edges = ncavity_edges + 1
                cavity_edges(1, ncavity_edges) = all_edges(1, i)
                cavity_edges(2, ncavity_edges) = all_edges(2, i)

                ! Find external neighbor (the triangle across this edge
                ! that is NOT in the cavity)
                t = edge_tri(i)
                j = edge_idx_in_tri(i)
                neighbor = mesh%triangles(t)%neighbors(j)

                ! Verify neighbor is not in cavity
                if (neighbor > 0) then
                    in_cavity = .false.
                    do edge_idx = 1, ncavity_triangles
                        if (cavity_triangles(edge_idx) == neighbor) then
                            in_cavity = .true.
                            exit
                        end if
                    end do
                    if (in_cavity) neighbor = 0
                end if

                external_neighbors(ncavity_edges) = neighbor
            end if
        end do
    end subroutine find_cavity_boundary

    subroutine remove_cavity_triangles(mesh, cavity_triangles, ncavity_triangles)
        !> Mark cavity triangles as invalid.
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: cavity_triangles(MAX_CAVITY_SIZE)
        integer, intent(in) :: ncavity_triangles

        integer :: i

        do i = 1, ncavity_triangles
            mesh%triangles(cavity_triangles(i))%valid = .false.
        end do
    end subroutine remove_cavity_triangles

    subroutine fill_cavity(mesh, point_idx, cavity_edges, external_neighbors,    &
                           ncavity_edges)
        !> Create new triangles with incremental adjacency updates.
        !
        !  Each new triangle has 3 edges:
        !    1. Boundary edge (v1-v2): neighbor is the external_neighbor
        !    2. Edge to next triangle in fan (v2-point): neighbor is next new tri
        !    3. Edge to previous triangle in fan (point-v1): neighbor is prev new tri
        !
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: point_idx
        integer, intent(in) :: cavity_edges(2, MAX_CAVITY_SIZE * 3)
        integer, intent(in) :: external_neighbors(MAX_CAVITY_SIZE * 3)
        integer, intent(in) :: ncavity_edges

        integer :: i, v1, v2, new_tri, ext_neighbor
        integer :: new_triangles(MAX_CAVITY_SIZE * 3)
        integer :: edge_v1(MAX_CAVITY_SIZE * 3)
        integer :: edge_v2(MAX_CAVITY_SIZE * 3)
        type(point_t) :: p, p1, p2
        integer :: j, k, other_v1, other_v2, edge_in_ext

        p = mesh%points(point_idx)

        ! First pass: create all new triangles
        do i = 1, ncavity_edges
            v1 = cavity_edges(1, i)
            v2 = cavity_edges(2, i)
            p1 = mesh%points(v1)
            p2 = mesh%points(v2)

            ! Ensure counter-clockwise orientation
            if (orientation(p1, p2, p) == ORIENTATION_CCW) then
                new_tri = add_triangle(mesh, v1, v2, point_idx)
                edge_v1(i) = v1
                edge_v2(i) = v2
            else
                new_tri = add_triangle(mesh, v2, v1, point_idx)
                edge_v1(i) = v2
                edge_v2(i) = v1
            end if
            new_triangles(i) = new_tri
        end do

        ! Second pass: set adjacency
        do i = 1, ncavity_edges
            new_tri = new_triangles(i)
            v1 = edge_v1(i)
            v2 = edge_v2(i)
            ext_neighbor = external_neighbors(i)

            ! Edge 1 (v1-v2): connects to external neighbor
            mesh%triangles(new_tri)%neighbors(1) = ext_neighbor

            ! Update external neighbor to point back to us
            if (ext_neighbor > 0 .and. ext_neighbor <= mesh%ntriangles) then
                call update_neighbor_link(mesh, ext_neighbor, v1, v2, new_tri)
            end if

            ! Edge 2 (v2-point): find which new triangle shares this edge
            do j = 1, ncavity_edges
                if (j == i) cycle
                other_v1 = edge_v1(j)
                other_v2 = edge_v2(j)
                ! Check if edge (v2, point_idx) matches edge (other_v1, point_idx)
                if (other_v1 == v2) then
                    mesh%triangles(new_tri)%neighbors(2) = new_triangles(j)
                    exit
                end if
            end do

            ! Edge 3 (point-v1): find which new triangle shares this edge
            do k = 1, ncavity_edges
                if (k == i) cycle
                other_v1 = edge_v1(k)
                other_v2 = edge_v2(k)
                ! Check if edge (point_idx, v1) matches edge (point_idx, other_v2)
                if (other_v2 == v1) then
                    mesh%triangles(new_tri)%neighbors(3) = new_triangles(k)
                    exit
                end if
            end do
        end do

        ! Update last_triangle for warm start
        if (ncavity_edges > 0) then
            last_triangle = new_triangles(1)
        end if
    end subroutine fill_cavity

    subroutine update_neighbor_link(mesh, tri_idx, v1, v2, new_neighbor)
        !> Update triangle tri_idx to point to new_neighbor across edge (v1,v2).
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: tri_idx, v1, v2, new_neighbor

        integer :: e, ev1, ev2

        do e = 1, 3
            ev1 = mesh%triangles(tri_idx)%vertices(e)
            ev2 = mesh%triangles(tri_idx)%vertices(mod(e, 3) + 1)

            ! Check if edge matches (in either direction)
            if ((ev1 == v1 .and. ev2 == v2) .or. (ev1 == v2 .and. ev2 == v1)) then
                mesh%triangles(tri_idx)%neighbors(e) = new_neighbor
                return
            end if
        end do
    end subroutine update_neighbor_link

end module bowyer_watson
