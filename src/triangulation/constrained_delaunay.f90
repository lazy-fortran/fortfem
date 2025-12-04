module constrained_delaunay
    !> Constrained Delaunay Triangulation
    !
    !  Algorithm based on the approach described in Shewchuks paper:
    !    1. Build unconstrained Delaunay triangulation (Bowyer-Watson)
    !    2. Insert constraint segments by flipping crossing edges
    !    3. After each flip, restore Delaunay property on affected triangles
    !    4. Remove exterior triangles using infecthull + plague algorithm
    !
    !  The key insight is that after flipping an edge to insert a constraint,
    !  we must recursively check and fix any edges that violate the Delaunay
    !  criterion (the incircle test). This ensures the final triangulation
    !  is a proper Constrained Delaunay Triangulation, not just a constrained
    !  triangulation.
    !
    use fortfem_kinds, only: dp
    use delaunay_types
    use geometric_predicates
    use bowyer_watson
    implicit none

    private
    public :: constrained_delaunay_triangulate

    ! Module-level storage for constraint segments (needed during fixup)
    integer, allocatable :: current_segments(:,:)

contains

    subroutine constrained_delaunay_triangulate(input_points,                 &
                                               constraint_segments, mesh)
        real(dp), intent(in) :: input_points(:,:)
        integer, intent(in) :: constraint_segments(:,:)
        type(mesh_t), intent(out) :: mesh

        integer :: i
        integer, allocatable :: adjusted_segments(:,:)
        integer :: offset

        call delaunay_triangulate(input_points, mesh)

        if (size(constraint_segments, 2) == 0) return

        ! Bowyer-Watson adds 3 super-triangle vertices at indices 1-3.
        ! The input points start at index 4. We need to offset the
        ! constraint segments accordingly.
        offset = 3
        allocate(adjusted_segments(2, size(constraint_segments, 2)))
        do i = 1, size(constraint_segments, 2)
            adjusted_segments(1, i) = constraint_segments(1, i) + offset
            adjusted_segments(2, i) = constraint_segments(2, i) + offset
        end do

        ! Store segments for use during Delaunay fixup
        if (allocated(current_segments)) deallocate(current_segments)
        allocate(current_segments, source=adjusted_segments)

        do i = 1, size(adjusted_segments, 2)
            call insert_segment(mesh, adjusted_segments(1, i),                &
                               adjusted_segments(2, i))
        end do

        call remove_exterior_triangles(mesh, adjusted_segments)

        deallocate(adjusted_segments)
        if (allocated(current_segments)) deallocate(current_segments)
    end subroutine constrained_delaunay_triangulate

    subroutine insert_segment(mesh, v1, v2)
        !> Insert a constraint segment by flipping edges that cross it.
        !  After each flip, restore Delaunay property on affected triangles.
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: v1, v2

        integer :: max_iter, iter
        logical :: segment_exists

        if (v1 == v2) return
        if (v1 < 1 .or. v1 > mesh%npoints) return
        if (v2 < 1 .or. v2 > mesh%npoints) return
        if (.not. mesh%points(v1)%valid) return
        if (.not. mesh%points(v2)%valid) return

        max_iter = 10 * mesh%ntriangles
        do iter = 1, max_iter
            segment_exists = edge_exists_in_mesh(mesh, v1, v2)
            if (segment_exists) exit
            if (.not. flip_one_crossing_edge(mesh, v1, v2)) exit
        end do
    end subroutine insert_segment

    logical function edge_exists_in_mesh(mesh, v1, v2) result(exists)
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: v1, v2

        integer :: t, a, b, c

        exists = .false.
        do t = 1, mesh%ntriangles
            if (.not. mesh%triangles(t)%valid) cycle
            a = mesh%triangles(t)%vertices(1)
            b = mesh%triangles(t)%vertices(2)
            c = mesh%triangles(t)%vertices(3)
            if (edges_match(a, b, v1, v2)) then
                exists = .true.
                return
            end if
            if (edges_match(b, c, v1, v2)) then
                exists = .true.
                return
            end if
            if (edges_match(c, a, v1, v2)) then
                exists = .true.
                return
            end if
        end do
    end function edge_exists_in_mesh

    logical function edges_match(a, b, v1, v2) result(match)
        integer, intent(in) :: a, b, v1, v2
        match = (a == v1 .and. b == v2) .or. (a == v2 .and. b == v1)
    end function edges_match

    logical function flip_one_crossing_edge(mesh, va, vb) result(flipped)
        !> Find and flip one edge that crosses segment va-vb.
        !  After the flip, restore Delaunay property on both sides.
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: va, vb

        integer :: t, i, e1, e2, t2, v_opp1, v_opp2
        integer :: new_t1, new_t2
        type(point_t) :: pa, pb, p1, p2, p3, p4
        real(dp) :: det1, det2

        flipped = .false.
        pa = mesh%points(va)
        pb = mesh%points(vb)

        do t = 1, mesh%ntriangles
            if (.not. mesh%triangles(t)%valid) cycle

            do i = 1, 3
                e1 = mesh%triangles(t)%vertices(i)
                e2 = mesh%triangles(t)%vertices(mod(i, 3) + 1)

                if (e1 == va .or. e1 == vb) cycle
                if (e2 == va .or. e2 == vb) cycle

                p1 = mesh%points(e1)
                p2 = mesh%points(e2)

                if (.not. segments_properly_intersect(pa, pb, p1, p2)) cycle

                t2 = find_adjacent_triangle(mesh, t, e1, e2)
                if (t2 == 0) cycle

                v_opp1 = get_opposite_vertex(mesh, t, e1, e2)
                v_opp2 = get_opposite_vertex(mesh, t2, e1, e2)
                if (v_opp1 == 0 .or. v_opp2 == 0) cycle
                if (v_opp1 == v_opp2) cycle

                p3 = mesh%points(v_opp1)
                p4 = mesh%points(v_opp2)

                ! Check if quadrilateral is convex (required for flip)
                det1 = signed_area(p3, p1, p4)
                det2 = signed_area(p3, p2, p4)
                if (det1 * det2 >= 0.0_dp) cycle

                det1 = signed_area(p4, p1, p3)
                det2 = signed_area(p4, p2, p3)
                if (det1 * det2 >= 0.0_dp) cycle

                call do_edge_flip_with_fixup(mesh, t, t2, e1, e2, v_opp1,     &
                                             v_opp2, va, vb, new_t1, new_t2)
                flipped = .true.
                return
            end do
        end do
    end function flip_one_crossing_edge

    subroutine do_edge_flip_with_fixup(mesh, t1, t2, e1, e2, v1, v2, seg_a,   &
                                       seg_b, new_t1, new_t2)
        !> Flip edge (e1,e2) to edge (v1,v2) and restore Delaunay property.
        !
        !  The edge (e1,e2) is shared by triangles t1 and t2.
        !  v1 is the opposite vertex in t1, v2 is the opposite vertex in t2.
        !  After flip, we have triangles (v1,e1,v2) and (v1,v2,e2).
        !  Then we call delaunay_fixup on edges opposite to seg_a and seg_b.
        !
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: t1, t2, e1, e2, v1, v2, seg_a, seg_b
        integer, intent(out) :: new_t1, new_t2

        type(point_t) :: p1, p2, pv1, pv2

        mesh%triangles(t1)%valid = .false.
        mesh%triangles(t2)%valid = .false.

        p1 = mesh%points(e1)
        p2 = mesh%points(e2)
        pv1 = mesh%points(v1)
        pv2 = mesh%points(v2)

        ! Create first new triangle: (v1, e1, v2) or (v1, v2, e1)
        if (signed_area(pv1, p1, pv2) > 0.0_dp) then
            new_t1 = add_triangle(mesh, v1, e1, v2)
        else
            new_t1 = add_triangle(mesh, v1, v2, e1)
        end if

        ! Create second new triangle: (v2, e2, v1) or (v2, v1, e2)
        if (signed_area(pv2, p2, pv1) > 0.0_dp) then
            new_t2 = add_triangle(mesh, v2, e2, v1)
        else
            new_t2 = add_triangle(mesh, v2, v1, e2)
        end if

        ! Restore Delaunay property on edges that dont cross the segment
        ! The constraint segment is (seg_a, seg_b).
        ! We need to fix edges on both sides of the new edge (v1, v2).
        call delaunay_fixup(mesh, new_t1, v1, e1, seg_a, seg_b, 0)
        call delaunay_fixup(mesh, new_t2, v2, e2, seg_a, seg_b, 0)
    end subroutine do_edge_flip_with_fixup

    recursive subroutine delaunay_fixup(mesh, tri, v_base, v_edge, seg_a,     &
                                        seg_b, depth)
        !> Restore Delaunay property for edges not crossing the constraint.
        !
        !  This is the key routine that makes CDT produce good triangles.
        !  After flipping an edge to insert a constraint segment, we check
        !  if any affected edges violate the Delaunay criterion (incircle test).
        !  If so, we flip them and recursively check the new edges.
        !
        !  Arguments:
        !    tri    - Triangle to check
        !    v_base - Vertex at the base of the edge we came from
        !    v_edge - Other vertex of the edge we came from
        !    seg_a, seg_b - The constraint segment endpoints
        !    depth  - Recursion depth (for safety limit)
        !
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: tri, v_base, v_edge, seg_a, seg_b, depth

        integer :: v_apex, adj_tri, v_far
        integer :: new_t1, new_t2
        type(point_t) :: p_base, p_edge, p_apex, p_far
        real(dp) :: incircle_result
        integer, parameter :: MAX_DEPTH = 100

        if (depth > MAX_DEPTH) return
        if (tri < 1 .or. tri > mesh%ntriangles) return
        if (.not. mesh%triangles(tri)%valid) return

        ! Find the apex vertex (the one thats not v_base or v_edge)
        v_apex = get_third_vertex(mesh, tri, v_base, v_edge)
        if (v_apex == 0) return

        ! Find adjacent triangle across edge (v_base, v_apex)
        adj_tri = find_adjacent_triangle(mesh, tri, v_base, v_apex)
        if (adj_tri == 0) return

        ! Find the far vertex in the adjacent triangle
        v_far = get_opposite_vertex(mesh, adj_tri, v_base, v_apex)
        if (v_far == 0) return
        if (v_far == v_edge) return  ! Degenerate case

        ! Dont flip constraint edges
        if (is_constraint_edge_static(v_base, v_apex)) return

        ! Get point coordinates
        p_base = mesh%points(v_base)
        p_edge = mesh%points(v_edge)
        p_apex = mesh%points(v_apex)
        p_far = mesh%points(v_far)

        ! Check if the triangles form a convex quadrilateral
        if (.not. is_convex_quad(p_edge, p_base, p_far, p_apex)) return

        ! Incircle test: is v_far inside the circumcircle of (v_edge, v_base,
        ! v_apex)?
        ! If positive, v_far is inside and we should flip.
        incircle_result = incircle_test(p_edge, p_base, p_apex, p_far)

        if (incircle_result > 0.0_dp) then
            ! Edge (v_base, v_apex) is not locally Delaunay - flip it
            call flip_edge_simple(mesh, tri, adj_tri, v_base, v_apex, v_edge, &
                                  v_far, new_t1, new_t2)

            ! Recursively fix the two new edges
            call delaunay_fixup(mesh, new_t1, v_base, v_edge, seg_a, seg_b,   &
                                depth + 1)
            call delaunay_fixup(mesh, new_t2, v_apex, v_edge, seg_a, seg_b,   &
                                depth + 1)
        end if
    end subroutine delaunay_fixup

    real(dp) function incircle_test(pa, pb, pc, pd) result(det)
        !> Incircle test: returns positive if pd is inside circumcircle of
        !  triangle (pa, pb, pc), negative if outside, zero if on circle.
        !
        !  The vertices pa, pb, pc must be in counterclockwise order.
        !
        type(point_t), intent(in) :: pa, pb, pc, pd

        real(dp) :: adx, ady, bdx, bdy, cdx, cdy
        real(dp) :: alift, blift, clift
        real(dp) :: bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady

        adx = pa%x - pd%x
        ady = pa%y - pd%y
        bdx = pb%x - pd%x
        bdy = pb%y - pd%y
        cdx = pc%x - pd%x
        cdy = pc%y - pd%y

        alift = adx * adx + ady * ady
        blift = bdx * bdx + bdy * bdy
        clift = cdx * cdx + cdy * cdy

        bdxcdy = bdx * cdy
        cdxbdy = cdx * bdy
        cdxady = cdx * ady
        adxcdy = adx * cdy
        adxbdy = adx * bdy
        bdxady = bdx * ady

        det = alift * (bdxcdy - cdxbdy)                                       &
            + blift * (cdxady - adxcdy)                                       &
            + clift * (adxbdy - bdxady)
    end function incircle_test

    logical function is_convex_quad(p1, p2, p3, p4) result(convex)
        !> Check if quadrilateral (p1, p2, p3, p4) is convex.
        !  Points are in order around the quadrilateral.
        type(point_t), intent(in) :: p1, p2, p3, p4

        real(dp) :: cross1, cross2

        ! For the quad to be convex, diagonals must intersect.
        ! Check that p1-p3 and p2-p4 properly intersect.
        cross1 = signed_area(p1, p3, p2)
        cross2 = signed_area(p1, p3, p4)
        if (cross1 * cross2 >= 0.0_dp) then
            convex = .false.
            return
        end if

        cross1 = signed_area(p2, p4, p1)
        cross2 = signed_area(p2, p4, p3)
        if (cross1 * cross2 >= 0.0_dp) then
            convex = .false.
            return
        end if

        convex = .true.
    end function is_convex_quad

    integer function get_third_vertex(mesh, tri, v1, v2) result(v3)
        !> Get the third vertex of a triangle given two vertices.
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: tri, v1, v2

        integer :: a, b, c

        v3 = 0
        if (tri < 1 .or. tri > mesh%ntriangles) return
        if (.not. mesh%triangles(tri)%valid) return

        a = mesh%triangles(tri)%vertices(1)
        b = mesh%triangles(tri)%vertices(2)
        c = mesh%triangles(tri)%vertices(3)

        if (a /= v1 .and. a /= v2) v3 = a
        if (b /= v1 .and. b /= v2) v3 = b
        if (c /= v1 .and. c /= v2) v3 = c
    end function get_third_vertex

    subroutine flip_edge_simple(mesh, t1, t2, e1, e2, v1, v2, new_t1, new_t2)
        !> Simple edge flip without constraint checking.
        !  Flip edge (e1,e2) to edge (v1,v2).
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: t1, t2, e1, e2, v1, v2
        integer, intent(out) :: new_t1, new_t2

        type(point_t) :: p1, p2, pv1, pv2

        mesh%triangles(t1)%valid = .false.
        mesh%triangles(t2)%valid = .false.

        p1 = mesh%points(e1)
        p2 = mesh%points(e2)
        pv1 = mesh%points(v1)
        pv2 = mesh%points(v2)

        if (signed_area(pv1, p1, pv2) > 0.0_dp) then
            new_t1 = add_triangle(mesh, v1, e1, v2)
        else
            new_t1 = add_triangle(mesh, v1, v2, e1)
        end if

        if (signed_area(pv2, p2, pv1) > 0.0_dp) then
            new_t2 = add_triangle(mesh, v2, e2, v1)
        else
            new_t2 = add_triangle(mesh, v2, v1, e2)
        end if
    end subroutine flip_edge_simple

    logical function is_constraint_edge_static(v1, v2) result(is_constraint)
        !> Check if edge (v1, v2) is a constraint edge using module storage.
        integer, intent(in) :: v1, v2

        integer :: i, a, b

        is_constraint = .false.
        if (.not. allocated(current_segments)) return

        do i = 1, size(current_segments, 2)
            a = current_segments(1, i)
            b = current_segments(2, i)
            if (edges_match(a, b, v1, v2)) then
                is_constraint = .true.
                return
            end if
        end do
    end function is_constraint_edge_static

    logical function segments_properly_intersect(p1, p2, q1, q2)              &
        result(intersect)
        type(point_t), intent(in) :: p1, p2, q1, q2

        real(dp) :: d1, d2, d3, d4

        intersect = .false.

        d1 = signed_area(p1, p2, q1)
        d2 = signed_area(p1, p2, q2)
        d3 = signed_area(q1, q2, p1)
        d4 = signed_area(q1, q2, p2)

        if (d1 * d2 >= 0.0_dp) return
        if (d3 * d4 >= 0.0_dp) return

        intersect = .true.
    end function segments_properly_intersect

    real(dp) function signed_area(pa, pb, pc)
        type(point_t), intent(in) :: pa, pb, pc
        signed_area = (pb%x - pa%x) * (pc%y - pa%y) -                         &
                      (pc%x - pa%x) * (pb%y - pa%y)
    end function signed_area

    integer function find_adjacent_triangle(mesh, tri, v1, v2) result(adj)
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: tri, v1, v2

        integer :: t, a, b, c

        adj = 0
        do t = 1, mesh%ntriangles
            if (t == tri) cycle
            if (.not. mesh%triangles(t)%valid) cycle
            a = mesh%triangles(t)%vertices(1)
            b = mesh%triangles(t)%vertices(2)
            c = mesh%triangles(t)%vertices(3)
            if (edges_match(a, b, v1, v2) .or. edges_match(b, c, v1, v2) .or. &
                edges_match(c, a, v1, v2)) then
                adj = t
                return
            end if
        end do
    end function find_adjacent_triangle

    integer function get_opposite_vertex(mesh, tri, v1, v2) result(opp)
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: tri, v1, v2

        integer :: a, b, c

        opp = 0
        a = mesh%triangles(tri)%vertices(1)
        b = mesh%triangles(tri)%vertices(2)
        c = mesh%triangles(tri)%vertices(3)
        if (a /= v1 .and. a /= v2) opp = a
        if (b /= v1 .and. b /= v2) opp = b
        if (c /= v1 .and. c /= v2) opp = c
    end function get_opposite_vertex

    subroutine remove_exterior_triangles(mesh, segments)
        !> Remove triangles outside the constrained boundary using
        !  infecthull + plague algorithm:
        !  1. infecthull: Mark triangles on convex hull not protected by
        !     constraint segments
        !  2. plague: Spread infection to neighbors across non-constraint edges
        !  3. Delete all infected triangles
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: segments(:,:)

        logical, allocatable :: infected(:)
        integer, allocatable :: virus_pool(:)
        integer :: pool_size, pool_head, pool_tail
        integer :: t, i, v1, v2, neighbor

        if (mesh%ntriangles == 0) return
        if (size(segments, 2) == 0) return

        allocate(infected(mesh%ntriangles))
        infected = .false.

        pool_size = mesh%ntriangles
        allocate(virus_pool(pool_size))
        pool_head = 1
        pool_tail = 0

        call infecthull(mesh, segments, infected, virus_pool, pool_tail)

        do while (pool_head <= pool_tail)
            t = virus_pool(pool_head)
            pool_head = pool_head + 1

            if (.not. mesh%triangles(t)%valid) cycle

            do i = 1, 3
                v1 = mesh%triangles(t)%vertices(i)
                v2 = mesh%triangles(t)%vertices(mod(i, 3) + 1)

                if (is_constraint_edge(v1, v2, segments)) cycle

                neighbor = find_adjacent_triangle(mesh, t, v1, v2)
                if (neighbor == 0) cycle
                if (infected(neighbor)) cycle
                if (.not. mesh%triangles(neighbor)%valid) cycle

                infected(neighbor) = .true.
                pool_tail = pool_tail + 1
                virus_pool(pool_tail) = neighbor
            end do
        end do

        do t = 1, mesh%ntriangles
            if (infected(t)) mesh%triangles(t)%valid = .false.
        end do

        deallocate(infected)
        deallocate(virus_pool)
    end subroutine remove_exterior_triangles

    subroutine infecthull(mesh, segments, infected, virus_pool, pool_tail)
        !> Mark triangles on convex hull that are not protected by segments.
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: segments(:,:)
        logical, intent(inout) :: infected(:)
        integer, intent(inout) :: virus_pool(:)
        integer, intent(inout) :: pool_tail

        integer :: t, i, v1, v2, neighbor
        logical :: is_hull_edge

        do t = 1, mesh%ntriangles
            if (.not. mesh%triangles(t)%valid) cycle
            if (infected(t)) cycle

            do i = 1, 3
                v1 = mesh%triangles(t)%vertices(i)
                v2 = mesh%triangles(t)%vertices(mod(i, 3) + 1)

                neighbor = find_adjacent_triangle(mesh, t, v1, v2)
                is_hull_edge = (neighbor == 0)

                if (is_hull_edge) then
                    if (.not. is_constraint_edge(v1, v2, segments)) then
                        infected(t) = .true.
                        pool_tail = pool_tail + 1
                        virus_pool(pool_tail) = t
                        exit
                    end if
                end if
            end do
        end do
    end subroutine infecthull

    logical function is_constraint_edge(v1, v2, segments) result(is_constraint)
        integer, intent(in) :: v1, v2
        integer, intent(in) :: segments(:,:)

        integer :: i, a, b

        is_constraint = .false.
        do i = 1, size(segments, 2)
            a = segments(1, i)
            b = segments(2, i)
            if (edges_match(a, b, v1, v2)) then
                is_constraint = .true.
                return
            end if
        end do
    end function is_constraint_edge

end module constrained_delaunay
