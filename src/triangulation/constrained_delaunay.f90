module constrained_delaunay
    use fortfem_kinds, only: dp
    use delaunay_types
    use geometric_predicates
    use bowyer_watson
    implicit none

    private
    public :: constrained_delaunay_triangulate
    public :: insert_constraint, recover_constraints
    public :: enforce_constraints
    public :: constraint_edge_exists

contains

    subroutine constrained_delaunay_triangulate(input_points,                 &
                                               constraint_segments, mesh)
        !> Constrained Delaunay triangulation using a Triangle-style
        !  algorithm: first build an unconstrained Delaunay triangulation
        !  with Bowyer-Watson, then insert each constraint segment by
        !  flipping intersecting edges until the segment appears.
        real(dp), intent(in) :: input_points(:,:)      ! (2, npoints)
        integer, intent(in) :: constraint_segments(:,:) ! (2, nsegments)
        type(mesh_t), intent(out) :: mesh

        integer :: i

        call delaunay_triangulate(input_points, mesh)

        if (size(constraint_segments, 2) <= 0) return

        do i = 1, size(constraint_segments, 2)
            call insert_constraint(mesh, constraint_segments(:, i))
        end do

        call enforce_constraints(mesh, constraint_segments)
    end subroutine constrained_delaunay_triangulate

    subroutine insert_constraint(mesh, constraint_edge)
        !> Insert a single constraint segment into the triangulation by
        !  iteratively flipping edges that intersect the segment until the
        !  segment itself becomes an edge of the triangulation.
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: constraint_edge(2)

        integer :: v1, v2
        integer :: max_iter, iter
        logical :: ok

        v1 = constraint_edge(1)
        v2 = constraint_edge(2)

        if (v1 == v2) return
        if (v1 < 1 .or. v1 > mesh%npoints) return
        if (v2 < 1 .or. v2 > mesh%npoints) return
        if (.not. mesh%points(v1)%valid) return
        if (.not. mesh%points(v2)%valid) return

        if (constraint_edge_exists(mesh, v1, v2)) then
            call add_constraint_edge(mesh, v1, v2)
            return
        end if

        max_iter = 1000
        do iter = 1, max_iter
            if (constraint_edge_exists(mesh, v1, v2)) exit
            ok = flip_one_intersecting_edge(mesh, v1, v2)
            if (.not. ok) exit
        end do

        if (constraint_edge_exists(mesh, v1, v2)) then
            call add_constraint_edge(mesh, v1, v2)
        end if
    end subroutine insert_constraint

    logical function constraint_edge_exists(mesh, v1, v2)
        !> Check if segment v1-v2 exists as an edge of any valid triangle.
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: v1, v2

        integer :: i, a, b, c

        constraint_edge_exists = .false.

        do i = 1, mesh%ntriangles
            if (.not. mesh%triangles(i)%valid) cycle

            a = mesh%triangles(i)%vertices(1)
            b = mesh%triangles(i)%vertices(2)
            c = mesh%triangles(i)%vertices(3)

            if (edge_matches(a, b, v1, v2)) then
                constraint_edge_exists = .true.
                return
            end if
            if (edge_matches(b, c, v1, v2)) then
                constraint_edge_exists = .true.
                return
            end if
            if (edge_matches(c, a, v1, v2)) then
                constraint_edge_exists = .true.
                return
            end if
        end do
    end function constraint_edge_exists

    logical function flip_one_intersecting_edge(mesh, v_start, v_end) result(done)
        !> Find a non-constrained edge that intersects the segment
        !  v_start-v_end and perform a flip of the two adjacent triangles.
        !  Returns true if a flip was performed.
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: v_start, v_end

        integer :: t, i, v1, v2, v3, t2, v4
        type(point_t) :: ps, pe, p1, p2

        done = .false.

        ps = mesh%points(v_start)
        pe = mesh%points(v_end)

        do t = 1, mesh%ntriangles
            if (.not. mesh%triangles(t)%valid) cycle

            do i = 1, 3
                v1 = mesh%triangles(t)%vertices(i)
                v2 = mesh%triangles(t)%vertices(mod(i, 3) + 1)

                if (v1 == v2) cycle

                if (v1 == v_start .or. v1 == v_end) cycle
                if (v2 == v_start .or. v2 == v_end) cycle

                if (edge_is_constrained(mesh, v1, v2)) cycle

                p1 = mesh%points(v1)
                p2 = mesh%points(v2)

                if (.not. segments_properly_intersect(ps, pe, p1, p2)) cycle

                t2 = find_adjacent_triangle(mesh, t, v1, v2)
                if (t2 == 0) cycle
                if (.not. mesh%triangles(t2)%valid) cycle

                v3 = third_vertex_of_triangle(mesh, t, v1, v2)
                v4 = third_vertex_of_triangle(mesh, t2, v1, v2)

                if (v3 == 0 .or. v4 == 0) cycle
                if (v3 == v4) cycle

                mesh%triangles(t)%valid = .false.
                mesh%triangles(t2)%valid = .false.

                call add_oriented_triangle(mesh, v3, v1, v4)
                call add_oriented_triangle(mesh, v4, v2, v3)

                done = .true.
                return
            end do
        end do
    end function flip_one_intersecting_edge

    logical function edge_matches(a, b, v1, v2) result(match)
        !> Check if edge (a,b) represents the same undirected edge as (v1,v2).
        integer, intent(in) :: a, b, v1, v2

        match = ((a == v1 .and. b == v2) .or. (a == v2 .and. b == v1))
    end function edge_matches

    logical function edge_is_constrained(mesh, v1, v2) result(is_constrained)
        !> Check if edge (v1,v2) has been marked as constrained in mesh%edges.
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: v1, v2

        integer :: e, a, b

        is_constrained = .false.

        do e = 1, mesh%nedges
            if (.not. mesh%edges(e)%valid) cycle
            if (.not. mesh%edges(e)%constrained) cycle

            a = mesh%edges(e)%endpoints(1)
            b = mesh%edges(e)%endpoints(2)

            if (edge_matches(a, b, v1, v2)) then
                is_constrained = .true.
                return
            end if
        end do
    end function edge_is_constrained

    integer function find_adjacent_triangle(mesh, tri_idx, v1, v2)           &
        result(adj_tri)
        !> Find the second triangle that shares edge (v1,v2) with tri_idx.
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: tri_idx, v1, v2

        integer :: t, a, b, c
        integer :: matches

        adj_tri = 0
        matches = 0

        do t = 1, mesh%ntriangles
            if (t == tri_idx) cycle
            if (.not. mesh%triangles(t)%valid) cycle

            a = mesh%triangles(t)%vertices(1)
            b = mesh%triangles(t)%vertices(2)
            c = mesh%triangles(t)%vertices(3)

            if (edge_matches(a, b, v1, v2) .or. edge_matches(b, c, v1, v2) &
                .or. edge_matches(c, a, v1, v2)) then
                adj_tri = t
                matches = matches + 1
            end if
        end do

        if (matches /= 1) then
            if (matches > 1) then
                write(*,'(A,4I6)') 'Warning: edge has ', matches,          &
                    ' adjacent triangles for vertices ', v1, v2
            end if
        end if
    end function find_adjacent_triangle

    integer function third_vertex_of_triangle(mesh, tri_idx, v1, v2)         &
        result(v3)
        !> Return the third vertex of triangle tri_idx that is not v1 or v2.
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: tri_idx, v1, v2

        integer :: a, b, c

        v3 = 0

        if (tri_idx < 1 .or. tri_idx > mesh%ntriangles) return
        if (.not. mesh%triangles(tri_idx)%valid) return

        a = mesh%triangles(tri_idx)%vertices(1)
        b = mesh%triangles(tri_idx)%vertices(2)
        c = mesh%triangles(tri_idx)%vertices(3)

        if (a /= v1 .and. a /= v2) then
            v3 = a
        else if (b /= v1 .and. b /= v2) then
            v3 = b
        else if (c /= v1 .and. c /= v2) then
            v3 = c
        end if
    end function third_vertex_of_triangle

    logical function segments_properly_intersect(p1, p2, q1, q2)             &
        result(intersects)
        !> Test if segments p1-p2 and q1-q2 intersect in their interiors
        !  (shared endpoints are not counted as intersections).
        type(point_t), intent(in) :: p1, p2, q1, q2

        integer :: o1, o2, o3, o4

        intersects = .false.

        o1 = orientation(p1, p2, q1)
        o2 = orientation(p1, p2, q2)
        o3 = orientation(q1, q2, p1)
        o4 = orientation(q1, q2, p2)

        if (o1 == ORIENTATION_COLLINEAR) return
        if (o2 == ORIENTATION_COLLINEAR) return
        if (o3 == ORIENTATION_COLLINEAR) return
        if (o4 == ORIENTATION_COLLINEAR) return

        if (o1 /= o2 .and. o3 /= o4) then
            intersects = .true.
        end if
    end function segments_properly_intersect

    subroutine add_oriented_triangle(mesh, v1, v2, v3)
        !> Add a triangle with vertices ordered to ensure positive area.
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: v1, v2, v3

        type(point_t) :: p1, p2, p3
        integer :: o
        integer :: tri_idx

        if (v1 == v2 .or. v2 == v3 .or. v3 == v1) return

        p1 = mesh%points(v1)
        p2 = mesh%points(v2)
        p3 = mesh%points(v3)

        o = orientation(p1, p2, p3)

        if (o == ORIENTATION_CCW) then
            tri_idx = add_triangle(mesh, v1, v2, v3)
        else if (o == ORIENTATION_CW) then
            tri_idx = add_triangle(mesh, v1, v3, v2)
        end if
    end subroutine add_oriented_triangle

    subroutine add_constraint_edge(mesh, v1, v2)
        !> Record a constraint edge in the mesh edge list.
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: v1, v2

        integer :: e

        e = add_edge(mesh, v1, v2, .true.)
    end subroutine add_constraint_edge

    subroutine recover_constraints(mesh, constraint_segments)
        !> Ensure all constraint segments appear as edges in the mesh.
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: constraint_segments(:,:)

        integer :: i

        do i = 1, size(constraint_segments, 2)
            if (.not. constraint_edge_exists(mesh, constraint_segments(1, i), &
                                             constraint_segments(2, i))) then
                call insert_constraint(mesh, constraint_segments(:, i))
            end if
        end do
    end subroutine recover_constraints

    subroutine enforce_constraints(mesh, constraint_segments)
        !> Final pass to reinsert any missing constraints (robustness loop).
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: constraint_segments(:,:)

        integer :: i, iter, max_iter
        logical :: all_ok

        max_iter = 10

        do iter = 1, max_iter
            all_ok = .true.
            do i = 1, size(constraint_segments, 2)
                if (.not. constraint_edge_exists(mesh,                     &
                                                constraint_segments(1, i), &
                                                constraint_segments(2, i))) then
                    call insert_constraint(mesh, constraint_segments(:, i))
                    all_ok = .false.
                end if
            end do
            if (all_ok) exit
        end do
    end subroutine enforce_constraints

end module constrained_delaunay
