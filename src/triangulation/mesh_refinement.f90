module mesh_refinement
    !> Delaunay refinement for quality meshing.
    !
    !  This module refines an existing constrained Delaunay triangulation by
    !  inserting Steiner points until all triangles satisfy a minimum angle
    !  constraint (and optionally a maximum area constraint).
    !
    !  Algorithm (Ruppert-style with Triangle-like processing):
    !    1. Split encroached segments: if any vertex lies inside the
    !       diametral circle of a constrained segment, split that segment at
    !       its midpoint and insert the new point with Bowyer-Watson.
    !    2. Split bad triangles: for any triangle whose minimum angle is
    !       below the target or whose area exceeds max_area, compute its
    !       circumcenter and insert it as a Steiner point. If the
    !       circumcenter would encroach a segment, split that segment first
    !       and retry the triangle later.
    !    3. Repeat until no encroached segments or bad triangles remain or
    !       a safety iteration limit is reached.
    !
    !  Key insight from Triangle library: process ALL bad triangles in each
    !  pass, not just one. When a circumcenter encroaches a segment, fix the
    !  segment but continue checking other triangles.
    !
    use fortfem_kinds, only: dp
    use delaunay_types, only: mesh_t, point_t, add_point
    use geometric_predicates, only: triangle_area, triangle_angles,           &
        circumcenter, geometric_tolerance, disable_robust_predicates
    use bowyer_watson, only: insert_point
    use point_location, only: build_adjacency, locate_point
    implicit none

    private
    public :: refine_delaunay

    real(dp), parameter :: PI = 3.14159265358979323846_dp

contains

    subroutine refine_delaunay(mesh, input_segments, min_angle_deg, max_area)
        !> Refine a constrained Delaunay mesh to improve triangle quality.
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: input_segments(:,:)
        real(dp), intent(in) :: min_angle_deg
        real(dp), intent(in) :: max_area

        integer, allocatable :: segments(:,:)
        integer :: nseg
        integer :: iter, max_iter
        logical :: changed

        if (mesh%ntriangles <= 0) return

        ! Disable robust predicates during refinement.
        ! Circumcenters are computed with float arithmetic, so the incircle
        ! test must also use float arithmetic to remain consistent.
        call disable_robust_predicates()

        nseg = size(input_segments, 2)
        if (nseg > 0) then
            allocate(segments(2, nseg))
            segments = input_segments
        end if

        ! Ensure adjacency is available for point insertion.
        call build_adjacency(mesh)

        max_iter = 1000
        do iter = 1, max_iter
            changed = .false.

            ! Phase 1: Fix all encroached segments first (like Triangle).
            call split_all_encroached_segments(mesh, segments, changed)

            ! Phase 2: Process bad triangles. Continue until no progress.
            call split_all_bad_triangles(mesh, segments, min_angle_deg,       &
                                         max_area, changed)

            if (.not. changed) exit
        end do

        if (allocated(segments)) deallocate(segments)
    end subroutine refine_delaunay

    subroutine split_all_encroached_segments(mesh, segments, changed)
        !> Split ALL encroached segments in one pass.
        type(mesh_t), intent(inout) :: mesh
        integer, allocatable, intent(inout) :: segments(:,:)
        logical, intent(inout) :: changed

        integer :: i, j, nseg
        integer :: v1, v2
        type(point_t) :: p1, p2, pv
        real(dp) :: cx, cy, dx, dy, r2, d2
        logical :: found_encroached
        logical :: seg_changed

        if (.not. allocated(segments)) return

        ! Keep splitting until no more encroached segments found.
        do
            found_encroached = .false.
            nseg = size(segments, 2)

            do i = 1, nseg
                v1 = segments(1, i)
                v2 = segments(2, i)
                if (v1 < 1 .or. v1 > mesh%npoints) cycle
                if (v2 < 1 .or. v2 > mesh%npoints) cycle
                if (.not. mesh%points(v1)%valid) cycle
                if (.not. mesh%points(v2)%valid) cycle

                p1 = mesh%points(v1)
                p2 = mesh%points(v2)

                cx = 0.5_dp * (p1%x + p2%x)
                cy = 0.5_dp * (p1%y + p2%y)
                dx = p1%x - p2%x
                dy = p1%y - p2%y
                r2 = 0.25_dp * (dx*dx + dy*dy) - geometric_tolerance
                if (r2 <= 0.0_dp) cycle

                do j = 1, mesh%npoints
                    if (.not. mesh%points(j)%valid) cycle
                    if (j == v1 .or. j == v2) cycle

                    pv = mesh%points(j)
                    dx = pv%x - cx
                    dy = pv%y - cy
                    d2 = dx*dx + dy*dy

                    if (d2 < r2) then
                        call split_segment_midpoint(mesh, segments, i,        &
                                                    seg_changed)
                        if (seg_changed) then
                            changed = .true.
                            found_encroached = .true.
                        end if
                        exit  ! Segment array changed, restart scan
                    end if
                end do

                if (found_encroached) exit  ! Restart from beginning
            end do

            if (.not. found_encroached) exit
        end do
    end subroutine split_all_encroached_segments

    subroutine split_all_bad_triangles(mesh, segments, min_angle_deg,         &
                                       max_area, changed)
        !> Process ALL bad triangles, not just one.
        !
        !  Uses Triangle-style off-center placement: when the circumcenter
        !  would fall outside the domain or too far from the triangle, we
        !  instead insert a point closer to the shortest edge. This
        !  guarantees termination without requiring segment splitting for
        !  boundary triangles.
        !
        type(mesh_t), intent(inout) :: mesh
        integer, allocatable, intent(inout) :: segments(:,:)
        real(dp), intent(in) :: min_angle_deg
        real(dp), intent(in) :: max_area
        logical, intent(inout) :: changed

        integer :: t, ntri_start
        integer :: v1, v2, v3
        type(point_t) :: pa, pb, pc
        type(point_t) :: steiner
        real(dp) :: area
        real(dp) :: angles(3)
        real(dp) :: tri_min
        integer :: enc_seg
        integer :: new_vertex
        integer :: seed_tri, loc_type
        logical :: seg_changed
        real(dp) :: offconstant

        if (mesh%ntriangles <= 0) return

        ! Compute off-center constant from minimum angle (Triangle formula).
        ! offconstant = 0.5 / tan(min_angle) ensures the off-center produces
        ! triangles with at least the target minimum angle.
        offconstant = 0.5_dp / tan(min_angle_deg * PI / 180.0_dp)

        ! Process triangles. Note: ntriangles may grow as we insert points,
        ! so we snapshot the count and only process existing triangles.
        ntri_start = mesh%ntriangles

        do t = 1, ntri_start
            if (.not. mesh%triangles(t)%valid) cycle

            v1 = mesh%triangles(t)%vertices(1)
            v2 = mesh%triangles(t)%vertices(2)
            v3 = mesh%triangles(t)%vertices(3)

            if (v1 < 1 .or. v1 > mesh%npoints) cycle
            if (v2 < 1 .or. v2 > mesh%npoints) cycle
            if (v3 < 1 .or. v3 > mesh%npoints) cycle
            if (.not. mesh%points(v1)%valid) cycle
            if (.not. mesh%points(v2)%valid) cycle
            if (.not. mesh%points(v3)%valid) cycle

            pa = mesh%points(v1)
            pb = mesh%points(v2)
            pc = mesh%points(v3)

            area = triangle_area(pa, pb, pc)
            if (area <= geometric_tolerance) cycle

            angles = triangle_angles(pa, pb, pc)
            tri_min = min(angles(1), min(angles(2), angles(3)))

            ! Skip good triangles
            if (tri_min >= min_angle_deg .and. area <= max_area) cycle

            ! Compute Steiner point using off-center placement (like Triangle).
            ! This guarantees the point is inside the domain for valid triangles.
            steiner = compute_offcenter(pa, pb, pc, offconstant)

            ! Check if Steiner point encroaches any segment
            enc_seg = find_encroached_segment_by_point(mesh, segments, steiner)
            if (enc_seg > 0) then
                ! Split the encroached segment instead of inserting Steiner point
                call split_segment_midpoint(mesh, segments, enc_seg,          &
                                            seg_changed)
                if (seg_changed) changed = .true.
                cycle
            end if

            ! Verify Steiner point is inside domain (should always be true with
            ! off-center, but check anyway for robustness)
            call locate_point(mesh, steiner%x, steiner%y, seed_tri, loc_type, 1)
            if (seed_tri == 0) then
                ! Off-center still outside - skip this triangle for now.
                ! It will be fixed when adjacent triangles are refined.
                cycle
            end if

            ! Insert the Steiner point
            new_vertex = add_point(mesh, steiner%x, steiner%y)
            call insert_point(mesh, new_vertex)

            ! Verify insertion succeeded
            if (.not. vertex_has_triangles(mesh, new_vertex)) then
                mesh%points(new_vertex)%valid = .false.
                cycle
            end if

            changed = .true.
        end do
    end subroutine split_all_bad_triangles

    function compute_offcenter(pa, pb, pc, offconstant) result(steiner)
        !> Compute off-center Steiner point for triangle (pa, pb, pc).
        !
        !  This implements the off-center placement from Triangle by Shewchuk.
        !  Instead of the true circumcenter (which may be far from the triangle
        !  for skinny triangles), we place the Steiner point closer to the
        !  shortest edge. The off-center is computed as:
        !
        !    midpoint of shortest edge + offconstant * perpendicular direction
        !
        !  where offconstant = 0.5 / tan(min_angle).
        !
        !  If the circumcenter is closer to the origin vertex than the off-center,
        !  we use the circumcenter instead.
        !
        type(point_t), intent(in) :: pa, pb, pc
        real(dp), intent(in) :: offconstant
        type(point_t) :: steiner

        real(dp) :: xdo, ydo, xao, yao
        real(dp) :: dodist, aodist, dadist
        real(dp) :: denominator
        real(dp) :: dx, dy, dxoff, dyoff

        ! Compute vectors from pa (origin) to pb (dest) and pc (apex)
        xdo = pb%x - pa%x
        ydo = pb%y - pa%y
        xao = pc%x - pa%x
        yao = pc%y - pa%y

        ! Compute squared edge lengths
        dodist = xdo * xdo + ydo * ydo  ! |pa-pb|^2
        aodist = xao * xao + yao * yao  ! |pa-pc|^2
        dadist = (pb%x - pc%x)**2 + (pb%y - pc%y)**2  ! |pb-pc|^2

        ! Compute circumcenter relative to pa
        denominator = 2.0_dp * (xdo * yao - xao * ydo)
        if (abs(denominator) < 1.0e-14_dp) then
            ! Degenerate triangle - return centroid
            steiner%x = (pa%x + pb%x + pc%x) / 3.0_dp
            steiner%y = (pa%y + pb%y + pc%y) / 3.0_dp
            steiner%id = 0
            return
        end if

        dx = (yao * dodist - ydo * aodist) / denominator
        dy = (xdo * aodist - xao * dodist) / denominator

        ! Find shortest edge and compute off-center if beneficial.
        ! The off-center is placed perpendicular to the shortest edge.
        if ((dodist < aodist) .and. (dodist < dadist)) then
            ! Edge pa-pb is shortest. Off-center perpendicular to it.
            dxoff = 0.5_dp * xdo - offconstant * ydo
            dyoff = 0.5_dp * ydo + offconstant * xdo
            ! Use off-center if closer to pa than circumcenter
            if (dxoff * dxoff + dyoff * dyoff < dx * dx + dy * dy) then
                dx = dxoff
                dy = dyoff
            end if
        else if (aodist < dadist) then
            ! Edge pa-pc is shortest. Off-center perpendicular to it.
            dxoff = 0.5_dp * xao + offconstant * yao
            dyoff = 0.5_dp * yao - offconstant * xao
            if (dxoff * dxoff + dyoff * dyoff < dx * dx + dy * dy) then
                dx = dxoff
                dy = dyoff
            end if
        else
            ! Edge pb-pc is shortest. Off-center perpendicular to it.
            dxoff = 0.5_dp * (pc%x - pb%x) - offconstant * (pc%y - pb%y)
            dyoff = 0.5_dp * (pc%y - pb%y) + offconstant * (pc%x - pb%x)
            ! Compare to circumcenter distance from pb
            if (dxoff * dxoff + dyoff * dyoff <                                &
                (dx - xdo) * (dx - xdo) + (dy - ydo) * (dy - ydo)) then
                dx = xdo + dxoff
                dy = ydo + dyoff
            end if
        end if

        steiner%x = pa%x + dx
        steiner%y = pa%y + dy
        steiner%id = 0
    end function compute_offcenter

    subroutine split_segment_midpoint(mesh, segments, seg_idx, changed)
        !> Split segment seg_idx at its midpoint and insert a new vertex.
        type(mesh_t), intent(inout) :: mesh
        integer, allocatable, intent(inout) :: segments(:,:)
        integer, intent(in) :: seg_idx
        logical, intent(out) :: changed

        integer :: v1, v2
        type(point_t) :: p1, p2, pmid
        integer :: new_vertex
        integer :: nseg, i, new_nseg
        integer, allocatable :: new_segments(:,:)

        changed = .false.
        if (.not. allocated(segments)) return

        nseg = size(segments, 2)
        if (seg_idx < 1 .or. seg_idx > nseg) return

        v1 = segments(1, seg_idx)
        v2 = segments(2, seg_idx)
        if (v1 < 1 .or. v1 > mesh%npoints) return
        if (v2 < 1 .or. v2 > mesh%npoints) return
        if (.not. mesh%points(v1)%valid) return
        if (.not. mesh%points(v2)%valid) return

        p1 = mesh%points(v1)
        p2 = mesh%points(v2)

        pmid%x = 0.5_dp * (p1%x + p2%x)
        pmid%y = 0.5_dp * (p1%y + p2%y)
        pmid%id = 0

        new_vertex = add_point(mesh, pmid%x, pmid%y)
        call insert_point(mesh, new_vertex)

        allocate(new_segments(2, nseg + 1))
        new_nseg = 0
        do i = 1, nseg
            if (i == seg_idx) cycle
            new_nseg = new_nseg + 1
            new_segments(1, new_nseg) = segments(1, i)
            new_segments(2, new_nseg) = segments(2, i)
        end do

        new_nseg = new_nseg + 1
        new_segments(1, new_nseg) = v1
        new_segments(2, new_nseg) = new_vertex

        new_nseg = new_nseg + 1
        new_segments(1, new_nseg) = new_vertex
        new_segments(2, new_nseg) = v2

        call move_alloc(new_segments, segments)

        changed = .true.
    end subroutine split_segment_midpoint

    logical function vertex_has_triangles(mesh, vertex_idx) result(has_tri)
        !> Check if a vertex is connected to any valid triangle.
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: vertex_idx

        integer :: t

        has_tri = .false.
        do t = 1, mesh%ntriangles
            if (.not. mesh%triangles(t)%valid) cycle
            if (any(mesh%triangles(t)%vertices == vertex_idx)) then
                has_tri = .true.
                return
            end if
        end do
    end function vertex_has_triangles

    integer function find_encroached_segment_by_point(mesh, segments, p)     &
        result(seg_idx)
        !> Find first segment whose diametral circle contains point p.
        type(mesh_t), intent(in) :: mesh
        integer, allocatable, intent(in) :: segments(:,:)
        type(point_t), intent(in) :: p

        integer :: i, nseg
        integer :: v1, v2
        type(point_t) :: p1, p2
        real(dp) :: cx, cy, dx, dy, r2, d2

        seg_idx = 0
        if (.not. allocated(segments)) return

        nseg = size(segments, 2)
        do i = 1, nseg
            v1 = segments(1, i)
            v2 = segments(2, i)
            if (v1 < 1 .or. v1 > mesh%npoints) cycle
            if (v2 < 1 .or. v2 > mesh%npoints) cycle
            if (.not. mesh%points(v1)%valid) cycle
            if (.not. mesh%points(v2)%valid) cycle

            p1 = mesh%points(v1)
            p2 = mesh%points(v2)

            cx = 0.5_dp * (p1%x + p2%x)
            cy = 0.5_dp * (p1%y + p2%y)
            dx = p1%x - p2%x
            dy = p1%y - p2%y
            r2 = 0.25_dp * (dx*dx + dy*dy) - geometric_tolerance
            if (r2 <= 0.0_dp) cycle

            dx = p%x - cx
            dy = p%y - cy
            d2 = dx*dx + dy*dy

            if (d2 < r2) then
                seg_idx = i
                return
            end if
        end do
    end function find_encroached_segment_by_point

end module mesh_refinement
