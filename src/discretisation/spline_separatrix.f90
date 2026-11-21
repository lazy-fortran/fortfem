module fortfem_spline_separatrix
    !! Block layouts aligned with separatrices through X-points inside the
    !! domain (doc/design/multipatch_spline_spaces.md, separatrix layouts).
    !!
    !! Every region bounded by separatrix branches is an O-grid about its
    !! O-point: a core around an inner polygon and a ring of blocks. Doublet
    !! lobes and the chain core have radial lines along rays from the
    !! O-point. An island is a lens between two branches from X_i to
    !! X_(i+1), and its O-grid follows the lens coordinates, so it stays
    !! inside a strongly curved island (layout_island_chain). The region
    !! between the outermost branches and the wall is a ring of blocks whose
    !! radial lines run to the wall, along the bisector of the outer sector
    !! at an X-point and along a ray from an O-point elsewhere. Every X-point
    !! is a vertex of valence 8. Both blocks at a branch piece carry
    !! LAYOUT_SEPARATRIX, and every block map is the Coons patch of its four
    !! sides. The topology does not depend on the samples, so moved curves
    !! keep the spline extraction (spline_space_update_geometry).
    use fortfem_kinds, only: dp
    use fortfem_spline_layout, only: block_layout_t, layout_curve_t, &
        layout_from_blocks, layout_curve_point, side_start, side_end, &
        LAYOUT_WALL, LAYOUT_SEPARATRIX
    use fortfem_spline_xpoint, only: sample_curve_t, closed_curve, open_curve, &
        ray_to_wall
    implicit none
    private

    type, public :: doublet_curves_t
        real(dp) :: xpoint(2) = 0.0_dp
        real(dp) :: opoint(2, 2) = 0.0_dp
        !> Loop k from X counter-clockwise about opoint(:, k) back to X. The
        !> first and last samples are replaced by X.
        type(sample_curve_t) :: lobe(2)
        !> Closed wall, counter-clockwise, not repeated at the end.
        real(dp), allocatable :: wall(:, :)
    end type doublet_curves_t

    type, public :: island_chain_curves_t
        !> Magnetic axis, inside the inner branches.
        real(dp) :: axis(2) = 0.0_dp
        !> X-points (2, m) counter-clockwise about the axis, and the island
        !> O-points (2, m), island i between X-points i and i + 1 (mod m).
        real(dp), allocatable :: xpoint(:, :), opoint(:, :)
        !> Branches from X-point i to X-point i + 1: inner(i) borders the
        !> core, outer(i) the outer region. End samples are replaced by the
        !> X-points.
        type(sample_curve_t), allocatable :: inner(:), outer(:)
        !> Closed wall, counter-clockwise, not repeated at the end.
        real(dp), allocatable :: wall(:, :)
    end type island_chain_curves_t

    type :: edge_t
        integer :: a = 0, b = 0, curve = 0
        real(dp) :: ua = 0.0_dp, ub = 0.0_dp
    end type edge_t

    type :: builder_t
        !! Vertices, curved edges and blocks of a layout under construction.
        !! Curve 1 is the wall. Sides on curves tagged LAYOUT_SEPARATRIX get
        !! that tag.
        real(dp), allocatable :: vertices(:, :)
        type(edge_t), allocatable :: edges(:)
        integer, allocatable :: corners(:, :)
        type(layout_curve_t), allocatable :: curves(:)
    end type builder_t

    public :: layout_doublet
    public :: layout_island_chain

contains

    subroutine layout_doublet(curves, core_fraction, layout, status, degree)
        !! Doublet: two lobes that touch at X, each a 5-block O-grid with
        !! outer corners X and the loop points at chord length 1/4, 1/2,
        !! 3/4, and a ring of 8 blocks to the wall, 18 blocks in all. Status
        !! 1: invalid input, 2: the lobes or the wall rays are not in
        !! counter-clockwise order about X, or a ray misses the wall, 3: a
        !! block is not counter-clockwise or the layout is not conforming.
        type(doublet_curves_t), intent(in) :: curves
        real(dp), intent(in) :: core_fraction
        type(block_layout_t), intent(out) :: layout
        integer, intent(out) :: status
        integer, intent(in), optional :: degree

        real(dp), parameter :: pi = acos(-1.0_dp)
        type(builder_t) :: b
        real(dp), allocatable :: points(:, :)
        real(dp) :: x(2), tangent(2), start(2), finish(2), turn
        real(dp) :: origin(2, 8), dir(2, 8), angle_start(2), angle_end(2)
        integer :: p, k, j, ix, outer(5, 2), ring(8)

        status = 1
        p = 3
        if (present(degree)) p = degree
        if (p < 1 .or. p > 7 .or. mod(p, 2) /= 1) return
        if (core_fraction <= 0.0_dp .or. core_fraction >= 1.0_dp) return
        if (.not. allocated(curves%wall)) return
        do k = 1, 2
            if (.not. allocated(curves%lobe(k)%points)) return
            if (size(curves%lobe(k)%points, 2) < 4) return
        end do
        allocate(b%curves(3), b%vertices(2, 0), b%edges(0), b%corners(4, 0))
        call closed_curve(curves%wall, p, LAYOUT_WALL, b%curves(1), status)
        if (status /= 0) return
        do k = 1, 2
            points = curves%lobe(k)%points
            points(:, 1) = curves%xpoint
            points(:, size(points, 2)) = curves%xpoint
            call open_curve(points, p, LAYOUT_SEPARATRIX, b%curves(1 + k), status)
            if (status /= 0) return
        end do
        status = 2
        do k = 1, 2
            if (loop_area(curves%lobe(k)%points) <= 0.0_dp) return
            call layout_curve_point(b%curves(1 + k), 0.0_dp, x, start)
            call layout_curve_point(b%curves(1 + k), 1.0_dp, x, finish)
            angle_start(k) = atan2(start(2), start(1))
            angle_end(k) = atan2(-finish(2), -finish(1))
        end do
        turn = sector(angle_start(1), angle_end(1)) + sector(angle_end(1), &
            angle_start(2)) + sector(angle_start(2), angle_end(2)) + &
            sector(angle_end(2), angle_start(1))
        if (abs(turn - 2.0_dp*pi) > 1.0e-8_dp) return

        call add_vertex(b, curves%xpoint, ix)
        do k = 1, 2
            outer([1, 5], k) = ix
            do j = 2, 4
                call layout_curve_point(b%curves(1 + k), 0.25_dp*(j - 1), x, tangent)
                call add_vertex(b, x, outer(j, k))
            end do
            do j = 1, 4
                call add_edge(b, outer(j, k), outer(j + 1, k), 1 + k, &
                    0.25_dp*(j - 1), 0.25_dp*j)
            end do
            call add_ogrid(b, curves%opoint(:, k), core_fraction, outer(:4, k))
        end do
        ! the ring runs counter-clockwise about X: the outer sector before
        ! lobe 1, lobe 1, the outer sector before lobe 2, lobe 2
        ring = [ix, outer(2:4, 1), ix, outer(2:4, 2)]
        origin = spread(curves%xpoint, 2, 8)
        dir(:, 1) = direction(angle_end(2) + 0.5_dp*sector(angle_end(2), &
            angle_start(1)))
        dir(:, 5) = direction(angle_end(1) + 0.5_dp*sector(angle_end(1), &
            angle_start(2)))
        do j = 2, 8
            if (j == 5) cycle
            k = merge(1, 2, j <= 4)
            origin(:, j) = curves%opoint(:, k)
            dir(:, j) = b%vertices(:, ring(j)) - origin(:, j)
        end do
        call add_wall_ring(b, ring, origin, dir, status)
        if (status /= 0) return
        call finish_layout(b, curves%xpoint, layout, status)
    end subroutine layout_doublet

    subroutine layout_island_chain(curves, core_fraction, layout, status, degree)
        !! Chain of m >= 2 islands about a magnetic axis. The core is an
        !! O-grid about the axis with outer corners X_1, M_1, ..., X_m, M_m
        !! (M_i the chord-length midpoint of inner branch i): a ring of 2m
        !! blocks around one quad (m = 2) or m quads about the axis as a
        !! vertex of valence m (m >= 3). Island i is a 5-block O-grid with
        !! outer corners X_i, the midpoint of outer branch i, X_(i+1) and M_i.
        !! A ring of 2m blocks runs from the outer branches to the wall:
        !! 19 blocks for m = 2, 10 m for m >= 3. Status 1: invalid input,
        !! 2: a branch loop or the wall rays are not counter-clockwise, or a
        !! ray misses the wall, 3: a block is not counter-clockwise or the
        !! layout is not conforming.
        type(island_chain_curves_t), intent(in) :: curves
        real(dp), intent(in) :: core_fraction
        type(block_layout_t), intent(out) :: layout
        integer, intent(out) :: status
        integer, intent(in), optional :: degree

        type(builder_t) :: b
        real(dp), allocatable :: points(:, :), origin(:, :), dir(:, :)
        real(dp) :: x(2), tangent(2), start(2), finish(2), d(2), axis(2), radius
        integer, allocatable :: xv(:), mi(:), me(:), core(:), inner(:), ring(:)
        integer :: p, m, i, in, j, jn, ia, io, c

        status = 1
        p = 3
        if (present(degree)) p = degree
        if (p < 1 .or. p > 7 .or. mod(p, 2) /= 1) return
        if (core_fraction <= 0.0_dp .or. core_fraction >= 1.0_dp) return
        if (.not. (allocated(curves%wall) .and. allocated(curves%xpoint) .and. &
            allocated(curves%opoint) .and. allocated(curves%inner) .and. &
            allocated(curves%outer))) return
        m = size(curves%xpoint, 2)
        if (m < 2 .or. any(shape(curves%opoint) /= [2, m]) .or. &
            size(curves%inner) /= m .or. size(curves%outer) /= m) return
        do i = 1, m
            if (.not. (allocated(curves%inner(i)%points) .and. &
                allocated(curves%outer(i)%points))) return
            if (size(curves%inner(i)%points, 2) < 3 .or. &
                size(curves%outer(i)%points, 2) < 3) return
        end do
        axis = curves%axis
        ! curve 1 + i: inner branch i, curve 1 + m + i: outer branch i
        allocate(b%curves(1 + 2*m), b%vertices(2, 0), b%edges(0), b%corners(4, 0))
        call closed_curve(curves%wall, p, LAYOUT_WALL, b%curves(1), status)
        if (status /= 0) return
        do i = 1, m
            in = modulo(i, m) + 1
            do j = 1, 2
                if (j == 1) then
                    points = curves%inner(i)%points
                else
                    points = curves%outer(i)%points
                end if
                points(:, 1) = curves%xpoint(:, i)
                points(:, size(points, 2)) = curves%xpoint(:, in)
                call open_curve(points, p, LAYOUT_SEPARATRIX, &
                    b%curves(1 + i + (j - 1)*m), status)
                if (status /= 0) return
            end do
        end do
        status = 2
        points = reshape([(curves%inner(i)%points, i = 1, m)], [2, &
            sum([(size(curves%inner(i)%points, 2), i = 1, m)])])
        if (loop_area(points) <= 0.0_dp) return
        do i = 1, m
            points = reshape([curves%outer(i)%points, &
                curves%inner(i)%points(:, size(curves%inner(i)%points, 2):1:-1)], &
                [2, size(curves%outer(i)%points, 2) + size(curves%inner(i)%points, 2)])
            if (loop_area(points) <= 0.0_dp) return
        end do

        allocate(xv(m), mi(m), me(m))
        do i = 1, m
            call add_vertex(b, curves%xpoint(:, i), xv(i))
        end do
        do i = 1, m
            in = modulo(i, m) + 1
            call layout_curve_point(b%curves(1 + i), 0.5_dp, x, tangent)
            call add_vertex(b, x, mi(i))
            call layout_curve_point(b%curves(1 + m + i), 0.5_dp, x, tangent)
            call add_vertex(b, x, me(i))
            call add_edge(b, xv(i), mi(i), 1 + i, 0.0_dp, 0.5_dp)
            call add_edge(b, mi(i), xv(in), 1 + i, 0.5_dp, 1.0_dp)
            call add_edge(b, xv(i), me(i), 1 + m + i, 0.0_dp, 0.5_dp)
            call add_edge(b, me(i), xv(in), 1 + m + i, 0.5_dp, 1.0_dp)
        end do

        ! core: ring of 2m blocks about the axis, then its inner polygon on
        ! the circle of core_fraction times the mean corner distance, which
        ! keeps the quads about the axis convex where the inner branches
        ! bulge towards the axis
        core = [(xv(i), mi(i), i = 1, m)]
        allocate(inner(2*m))
        radius = core_fraction*sum(norm2(b%vertices(:, core) - spread(axis, 2, 2*m), &
            dim=1))/(2*m)
        do j = 1, 2*m
            d = b%vertices(:, core(j)) - axis
            call add_vertex(b, axis + radius*d/norm2(d), inner(j))
        end do
        do j = 1, 2*m
            jn = modulo(j, 2*m) + 1
            call add_block(b, [core(j), core(jn), inner(jn), inner(j)])
        end do
        if (m == 2) then
            call add_block(b, inner)
        else
            call add_vertex(b, axis, ia)
            do i = 1, m
                c = 2*i - 1
                call add_block(b, [ia, inner(c), inner(c + 1), &
                    inner(modulo(c + 1, 2*m) + 1)])
            end do
        end if
        do i = 1, m
            in = modulo(i, m) + 1
            call add_lens(b, xv(i), xv(in), me(i), mi(i), 1 + i, 1 + m + i, &
                core_fraction, p, status)
            if (status /= 0) return
        end do

        ! outer ring: X_1, E_1 midpoint, X_2, ... counter-clockwise
        ring = [(xv(i), me(i), i = 1, m)]
        allocate(origin(2, 2*m), dir(2, 2*m))
        do i = 1, m
            io = merge(m, i - 1, i == 1)
            call layout_curve_point(b%curves(1 + m + io), 1.0_dp, x, finish)
            call layout_curve_point(b%curves(1 + m + i), 0.0_dp, x, start)
            d = start/norm2(start) - finish/norm2(finish)
            if (norm2(d) < 1.0e-8_dp) then
                d = curves%xpoint(:, i) - axis
            else if (dot_product(d, curves%xpoint(:, i) - axis) < 0.0_dp) then
                d = -d
            end if
            origin(:, 2*i - 1) = curves%xpoint(:, i)
            dir(:, 2*i - 1) = d
            origin(:, 2*i) = axis
            dir(:, 2*i) = b%vertices(:, me(i)) - axis
        end do
        call add_wall_ring(b, ring, origin, dir, status)
        if (status /= 0) return
        call finish_layout(b, axis, layout, status)
    end subroutine layout_island_chain

    subroutine add_ogrid(b, center, core_fraction, outer)
        !! Core quad and ring of four blocks inside the counter-clockwise
        !! outer corners, the inner corners on the rays from center.
        type(builder_t), intent(inout) :: b
        real(dp), intent(in) :: center(2), core_fraction
        integer, intent(in) :: outer(4)

        integer :: inner(4), j, jn

        do j = 1, 4
            call add_vertex(b, center + core_fraction*(b%vertices(:, outer(j)) - &
                center), inner(j))
        end do
        call add_block(b, inner)
        do j = 1, 4
            jn = modulo(j, 4) + 1
            call add_block(b, [outer(j), outer(jn), inner(jn), inner(j)])
        end do
    end subroutine add_ogrid

    subroutine add_lens(b, xa, xb, me, mi, ci, ce, core_fraction, p, status)
        !! O-grid of the lens between the inner branch ci and the outer branch
        !! ce from xa to xb, in the coordinates P(u, v) = (I(u) + E(u))/2 +
        !! v (E(u) - I(u))/2 at the common chord-length parameter u. The
        !! core corners are P(d, 0), P(1/2, c), P(1 - d, 0), P(1/2, -c),
        !! c = core_fraction, d = (1 - c)/2, and all interior sides follow P,
        !! so the blocks stay inside a strongly curved island.
        type(builder_t), intent(inout) :: b
        integer, intent(in) :: xa, xb, me, mi, ci, ce, p
        real(dp), intent(in) :: core_fraction
        integer, intent(out) :: status

        integer, parameter :: n = 21
        real(dp) :: c, d, path(2, 4), y(2)
        integer :: corner(4), side(4), j, k

        c = core_fraction
        d = 0.5_dp*(1.0_dp - c)
        ! (u, v) of the core corners: towards xa, me, xb, mi
        path = reshape([d, 0.0_dp, 0.5_dp, c, 1.0_dp - d, 0.0_dp, 0.5_dp, -c], [2, 4])
        do k = 1, 4
            y = lens_point(path(1, k), path(2, k))
            call add_vertex(b, y, corner(k))
        end do
        do k = 1, 4
            j = modulo(k, 4) + 1
            call add_lens_curve(path(1, k), path(2, k), path(1, j), path(2, j), &
                side(k))
            if (status /= 0) return
            call add_edge(b, corner(k), corner(j), side(k), 0.0_dp, 1.0_dp)
        end do
        ! every interior side is its own curve: a curve is parametrised by
        ! chord length, not by u
        call add_lens_curve(0.0_dp, 0.0_dp, d, 0.0_dp, j)
        if (status /= 0) return
        call add_edge(b, xa, corner(1), j, 0.0_dp, 1.0_dp)
        call add_lens_curve(1.0_dp, 0.0_dp, 1.0_dp - d, 0.0_dp, j)
        if (status /= 0) return
        call add_edge(b, xb, corner(3), j, 0.0_dp, 1.0_dp)
        call add_lens_curve(0.5_dp, 1.0_dp, 0.5_dp, c, j)
        if (status /= 0) return
        call add_edge(b, me, corner(2), j, 0.0_dp, 1.0_dp)
        call add_lens_curve(0.5_dp, -1.0_dp, 0.5_dp, -c, j)
        if (status /= 0) return
        call add_edge(b, mi, corner(4), j, 0.0_dp, 1.0_dp)
        call add_block(b, corner)
        call add_block(b, [xa, me, corner(2), corner(1)])
        call add_block(b, [me, xb, corner(3), corner(2)])
        call add_block(b, [xb, mi, corner(4), corner(3)])
        call add_block(b, [mi, xa, corner(1), corner(4)])

    contains

        function lens_point(u, v) result(y)
            real(dp), intent(in) :: u, v
            real(dp) :: y(2), yi(2), ye(2), t(2)

            call layout_curve_point(b%curves(ci), u, yi, t)
            call layout_curve_point(b%curves(ce), u, ye, t)
            y = 0.5_dp*(yi + ye) + 0.5_dp*v*(ye - yi)
        end function lens_point

        subroutine add_lens_curve(u0, v0, u1, v1, index)
            !! Interior curve through P on the straight (u, v) path.
            real(dp), intent(in) :: u0, v0, u1, v1
            integer, intent(out) :: index
            type(layout_curve_t) :: curve
            real(dp) :: samples(2, n), t
            integer :: q

            do q = 1, n
                t = real(q - 1, dp)/(n - 1)
                samples(:, q) = lens_point(u0 + t*(u1 - u0), v0 + t*(v1 - v0))
            end do
            call open_curve(samples, p, 0, curve, status)
            if (status /= 0) return
            call add_curve(b, curve, index)
        end subroutine add_lens_curve

    end subroutine add_lens

    subroutine add_curve(b, curve, index)
        type(builder_t), intent(inout) :: b
        type(layout_curve_t), intent(in) :: curve
        integer, intent(out) :: index
        type(layout_curve_t), allocatable :: grown(:)

        allocate(grown(size(b%curves) + 1))
        grown(:size(b%curves)) = b%curves
        grown(size(grown)) = curve
        call move_alloc(grown, b%curves)
        index = size(b%curves)
    end subroutine add_curve

    subroutine add_wall_ring(b, ring, origin, dir, status)
        !! Blocks between the closed counter-clockwise vertex cycle ring and
        !! the wall, with radial lines along the rays origin + t dir.
        type(builder_t), intent(inout) :: b
        integer, intent(in) :: ring(:)
        real(dp), intent(in) :: origin(:, :), dir(:, :)
        integer, intent(out) :: status

        real(dp) :: uwall(size(ring)), x(2), tangent(2), turn
        integer :: w(size(ring)), n, j, jn

        n = size(ring)
        do j = 1, n
            call ray_to_wall(b%curves(1), origin(:, j), origin(:, j) + dir(:, j), &
                uwall(j), status)
            if (status /= 0) return
            call layout_curve_point(b%curves(1), uwall(j), x, tangent)
            call add_vertex(b, x, w(j))
        end do
        status = 2
        turn = 0.0_dp
        do j = 1, n
            jn = modulo(j, n) + 1
            if (modulo(uwall(jn) - uwall(j), 1.0_dp) <= 0.0_dp) return
            turn = turn + modulo(uwall(jn) - uwall(j), 1.0_dp)
        end do
        if (abs(turn - 1.0_dp) > 1.0e-8_dp) return
        do j = 1, n
            jn = modulo(j, n) + 1
            call add_edge(b, w(j), w(jn), 1, uwall(j), uwall(j) + &
                modulo(uwall(jn) - uwall(j), 1.0_dp))
            call add_block(b, [ring(j), w(j), w(jn), ring(jn)])
        end do
        status = 0
    end subroutine add_wall_ring

    subroutine finish_layout(b, center, layout, status)
        !! Connect the blocks, attach the curve windows to their sides and
        !! tag the branch pieces.
        type(builder_t), intent(in) :: b
        real(dp), intent(in) :: center(2)
        type(block_layout_t), intent(out) :: layout
        integer, intent(out) :: status

        integer, allocatable :: tags(:, :)
        integer :: j, s, nb, ierr

        nb = size(b%corners, 2)
        allocate(tags(4, nb))
        tags = LAYOUT_WALL
        status = 3
        call layout_from_blocks(b%vertices, b%corners, tags, layout, ierr, &
            oriented=.true.)
        if (ierr /= 0) return
        layout%curves = b%curves
        layout%coons = .true.
        layout%center = center
        do j = 1, nb
            do s = 1, 4
                call side_window(b, b%corners(side_start(s), j), &
                    b%corners(side_end(s), j), layout%side_curve(s, j), &
                    layout%side_u0(s, j), layout%side_du(s, j))
                if (layout%side_curve(s, j) == 0) cycle
                if (layout%curves(layout%side_curve(s, j))%tag == LAYOUT_SEPARATRIX) &
                    layout%side_tag(s, j) = LAYOUT_SEPARATRIX
            end do
        end do
        status = 0
    end subroutine finish_layout

    subroutine add_vertex(b, y, index)
        type(builder_t), intent(inout) :: b
        real(dp), intent(in) :: y(2)
        integer, intent(out) :: index

        b%vertices = reshape([b%vertices, y], [2, size(b%vertices, 2) + 1])
        index = size(b%vertices, 2)
    end subroutine add_vertex

    subroutine add_edge(b, a, e, c, ua, ub)
        type(builder_t), intent(inout) :: b
        integer, intent(in) :: a, e, c
        real(dp), intent(in) :: ua, ub
        type(edge_t), allocatable :: grown(:)

        allocate(grown(size(b%edges) + 1))
        grown(:size(b%edges)) = b%edges
        grown(size(grown)) = edge_t(a, e, c, ua, ub)
        call move_alloc(grown, b%edges)
    end subroutine add_edge

    subroutine add_block(b, c)
        !! Corners in either orientation, stored counter-clockwise. The
        !! orientation is that of the boundary through the curved sides,
        !! since a thin curved block can have a clockwise corner polygon.
        type(builder_t), intent(inout) :: b
        integer, intent(in) :: c(4)

        integer, parameter :: n = 16
        real(dp) :: boundary(2, 4*n), t(2)
        real(dp) :: u0, du
        integer :: j, q, curve

        do j = 1, 4
            call side_window(b, c(j), c(modulo(j, 4) + 1), curve, u0, du)
            do q = 1, n
                if (curve == 0) then
                    boundary(:, n*(j - 1) + q) = b%vertices(:, c(j)) + &
                        real(q - 1, dp)/n*(b%vertices(:, c(modulo(j, 4) + 1)) - &
                        b%vertices(:, c(j)))
                else
                    call layout_curve_point(b%curves(curve), b%curves(curve)%offset + &
                        u0 + du*real(q - 1, dp)/n, boundary(:, n*(j - 1) + q), t)
                end if
            end do
        end do
        if (loop_area(boundary) >= 0.0_dp) then
            b%corners = reshape([b%corners, c], [4, size(b%corners, 2) + 1])
        else
            b%corners = reshape([b%corners, c(1), c(4), c(3), c(2)], &
                [4, size(b%corners, 2) + 1])
        end if
    end subroutine add_block

    subroutine side_window(b, a, e, c, u0, du)
        type(builder_t), intent(in) :: b
        integer, intent(in) :: a, e
        integer, intent(out) :: c
        real(dp), intent(out) :: u0, du
        integer :: k

        c = 0
        u0 = 0.0_dp
        du = 0.0_dp
        do k = 1, size(b%edges)
            if (b%edges(k)%a == a .and. b%edges(k)%b == e) then
                c = b%edges(k)%curve
                u0 = b%edges(k)%ua
                du = b%edges(k)%ub - b%edges(k)%ua
                return
            else if (b%edges(k)%a == e .and. b%edges(k)%b == a) then
                c = b%edges(k)%curve
                u0 = b%edges(k)%ub
                du = b%edges(k)%ua - b%edges(k)%ub
                return
            end if
        end do
    end subroutine side_window

    pure real(dp) function sector(from, to)
        !! Counter-clockwise angle from one direction to another, in [0, 2 pi).
        real(dp), intent(in) :: from, to

        sector = modulo(to - from, 2.0_dp*acos(-1.0_dp))
    end function sector

    pure function direction(angle) result(d)
        real(dp), intent(in) :: angle
        real(dp) :: d(2)

        d = [cos(angle), sin(angle)]
    end function direction

    pure real(dp) function loop_area(p) result(area)
        !! Signed area of the closed polygon through the columns of p.
        real(dp), intent(in) :: p(:, :)
        integer :: i, j

        area = 0.0_dp
        do i = 1, size(p, 2)
            j = modulo(i, size(p, 2)) + 1
            area = area + p(1, i)*p(2, j) - p(1, j)*p(2, i)
        end do
        area = 0.5_dp*area
    end function loop_area

end module fortfem_spline_separatrix
