module fortfem_spline_xpoint
    !! Block layouts with X-points on the plasma boundary
    !! (doc/design/multipatch_spline_spaces.md, X-point layouts).
    !!
    !! The core boundary Gamma_p has four corners K1 (bottom), K2 (outboard),
    !! K3 (top), K4 (inboard), counter-clockwise. A corner with legs is an
    !! X-point: its two legs run to the wall and bound a private flux region
    !! (three quads around a centre vertex). A smooth corner has a straight
    !! line along the ray from the centre to the wall. The core is a 5-block
    !! O-grid, and the scrape-off layer of each core segment is split into
    !! an inner and an outer layer at the mid points of its two radial
    !! sides. Every block map is the Coons patch of its four sides.
    use fortfem_kinds, only: dp
    use fortfem_spline_layout, only: block_layout_t, layout_curve_t, &
        layout_from_blocks, layout_curve_point, side_start, side_end, &
        LAYOUT_WALL, LAYOUT_PLASMA_BOUNDARY, LAYOUT_SEPARATRIX
    use fortfem_spline_extraction_1d, only: build_periodic_extraction_1d, &
        local_bsplines
    implicit none
    private

    type, public :: sample_curve_t
        real(dp), allocatable :: points(:, :) ! (2, n)
    end type sample_curve_t

    type, public :: xpoint_curves_t
        !> Gamma_p from corner q to corner q + 1, counter-clockwise.
        type(sample_curve_t) :: segment(4)
        !> Legs from corner q to the wall; present at X-points. leg_after
        !> borders the scrape-off layer of segment q, leg_before that of
        !> segment q - 1.
        type(sample_curve_t) :: leg_after(4), leg_before(4)
        !> Closed wall, counter-clockwise, not repeated at the end.
        real(dp), allocatable :: wall(:, :)
    end type xpoint_curves_t

    type :: edge_t
        integer :: a = 0, b = 0, curve = 0
        real(dp) :: ua = 0.0_dp, ub = 0.0_dp
    end type edge_t

    public :: layout_xpoint

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

contains

    subroutine layout_xpoint(center, core_fraction, curves, layout, status, degree)
        !! Layout from sampled Gamma_p segments, legs and wall. The same
        !! corners with legs give the same topology, so a new shape keeps the
        !! spline extraction (spline_space_update_geometry). Status 1:
        !! invalid input, 2: a leg or ray does not reach the wall, 3: a block
        !! is not counter-clockwise or the layout is not conforming.
        real(dp), intent(in) :: center(2), core_fraction
        type(xpoint_curves_t), intent(in) :: curves
        type(block_layout_t), intent(out) :: layout
        integer, intent(out) :: status
        integer, intent(in), optional :: degree

        type(layout_curve_t), allocatable :: curve(:)
        type(edge_t), allocatable :: edges(:)
        real(dp), allocatable :: vertices(:, :), points(:, :)
        integer, allocatable :: corners(:, :), tags(:, :)
        real(dp) :: uwall(4, 2), tangent(2), x(2), d, mid
        integer :: p, q, qn, nc, iw, ma(4), mb(4), ea(4), eb(4)
        integer :: k(4), inner(4), mw(4), cc(4), la(4), lb(4), nb, s, i, j
        logical :: xpt(4)

        status = 1
        p = 3
        if (present(degree)) p = degree
        if (p < 1 .or. p > 7 .or. mod(p, 2) /= 1) return
        if (core_fraction <= 0.0_dp .or. core_fraction >= 1.0_dp) return
        if (.not. allocated(curves%wall)) return
        do q = 1, 4
            if (.not. allocated(curves%segment(q)%points)) return
            xpt(q) = allocated(curves%leg_after(q)%points)
            if (xpt(q) .neqv. allocated(curves%leg_before(q)%points)) return
        end do
        allocate(curve(5 + 2*count(xpt)))
        call closed_curve(curves%wall, p, LAYOUT_WALL, curve(1), status)
        if (status /= 0) return
        iw = 1
        nc = 1
        allocate(vertices(2, 0))
        do q = 1, 4
            call add_vertex(curves%segment(q)%points(:, 1), k(q))
        end do
        do q = 1, 4
            qn = modulo(q, 4) + 1
            points = curves%segment(q)%points
            points(:, 1) = vertices(:, k(q))
            points(:, size(points, 2)) = vertices(:, k(qn))
            nc = nc + 1
            call open_curve(points, p, LAYOUT_PLASMA_BOUNDARY, curve(nc), status)
            if (status /= 0) return
        end do
        do q = 1, 4
            call add_vertex(center + core_fraction*(vertices(:, k(q)) - center), &
                inner(q))
        end do
        status = 2
        do q = 1, 4
            if (xpt(q)) then
                call leg(curves%leg_after(q)%points, uwall(q, 1), la(q))
                if (status /= 0) return
                call leg(curves%leg_before(q)%points, uwall(q, 2), lb(q))
                if (status /= 0) return
                call layout_curve_point(curve(la(q)), 1.0_dp, x, tangent)
                call add_vertex(x, ea(q))
                call layout_curve_point(curve(lb(q)), 1.0_dp, x, tangent)
                call add_vertex(x, eb(q))
                call layout_curve_point(curve(la(q)), 0.5_dp, x, tangent)
                call add_vertex(x, ma(q))
                call layout_curve_point(curve(lb(q)), 0.5_dp, x, tangent)
                call add_vertex(x, mb(q))
            else
                la(q) = 0
                lb(q) = 0
                call ray_to_wall(curve(iw), center, vertices(:, k(q)), uwall(q, 1), &
                    status)
                if (status /= 0) return
                uwall(q, 2) = uwall(q, 1)
                call layout_curve_point(curve(iw), uwall(q, 1), x, tangent)
                call add_vertex(x, ea(q))
                eb(q) = ea(q)
                call add_vertex(0.5_dp*(vertices(:, k(q)) + x), ma(q))
                mb(q) = ma(q)
            end if
        end do
        do q = 1, 4
            mw(q) = 0
            cc(q) = 0
            if (.not. xpt(q)) cycle
            d = modulo(uwall(q, 1) - uwall(q, 2), 1.0_dp)
            mid = uwall(q, 1) - 0.5_dp*d
            call layout_curve_point(curve(iw), mid, x, tangent)
            call add_vertex(x, mw(q))
            call add_vertex((vertices(:, ma(q)) + vertices(:, mb(q)) + x)/3.0_dp, cc(q))
        end do
        allocate(edges(0))
        do q = 1, 4
            qn = modulo(q, 4) + 1
            call add_edge(k(q), k(qn), 1 + q, 0.0_dp, 1.0_dp)
            call add_edge(k(q), ma(q), la(q), 0.0_dp, 0.5_dp)
            call add_edge(ma(q), ea(q), la(q), 0.5_dp, 1.0_dp)
            if (xpt(q)) then
                call add_edge(k(q), mb(q), lb(q), 0.0_dp, 0.5_dp)
                call add_edge(mb(q), eb(q), lb(q), 0.5_dp, 1.0_dp)
                d = modulo(uwall(q, 1) - uwall(q, 2), 1.0_dp)
                call add_edge(ea(q), mw(q), iw, uwall(q, 1), uwall(q, 1) - 0.5_dp*d)
                call add_edge(mw(q), eb(q), iw, uwall(q, 1) - 0.5_dp*d, &
                    uwall(q, 1) - d)
            end if
            d = modulo(uwall(qn, 2) - uwall(q, 1), 1.0_dp)
            call add_edge(ea(q), eb(qn), iw, uwall(q, 1), uwall(q, 1) + d)
        end do
        nb = 5 + 8 + 3*count(xpt)
        allocate(corners(4, nb))
        corners(:, 1) = inner
        nb = 1
        do q = 1, 4
            qn = modulo(q, 4) + 1
            call add_block([k(q), k(qn), inner(qn), inner(q)])
        end do
        do q = 1, 4
            qn = modulo(q, 4) + 1
            call add_block([ma(q), mb(qn), k(qn), k(q)])
            call add_block([ea(q), eb(qn), mb(qn), ma(q)])
        end do
        do q = 1, 4
            if (.not. xpt(q)) cycle
            call add_block([k(q), ma(q), cc(q), mb(q)])
            call add_block([ma(q), ea(q), mw(q), cc(q)])
            call add_block([cc(q), mw(q), eb(q), mb(q)])
        end do
        allocate(tags(4, nb))
        tags = LAYOUT_WALL
        status = 3
        call layout_from_blocks(vertices, corners(:, :nb), tags, layout, i)
        if (i /= 0) return
        layout%curves = curve(:nc)
        layout%coons(2:) = .true.
        layout%center = center
        do j = 1, nb
            do s = 1, 4
                call side_window(corners(side_start(s), j), corners(side_end(s), j), &
                    layout%side_curve(s, j), layout%side_u0(s, j), &
                    layout%side_du(s, j))
                if (layout%side_curve(s, j) == 0) cycle
                if (layout%curves(layout%side_curve(s, j))%tag == &
                    LAYOUT_PLASMA_BOUNDARY .and. j <= 5) &
                    layout%side_tag(s, j) = LAYOUT_PLASMA_BOUNDARY
                if (layout%curves(layout%side_curve(s, j))%tag == LAYOUT_SEPARATRIX &
                    .and. j >= 6 .and. j <= 13) &
                    layout%side_tag(s, j) = LAYOUT_SEPARATRIX
            end do
        end do
        status = 0

    contains

        subroutine add_vertex(x, index)
            real(dp), intent(in) :: x(2)
            integer, intent(out) :: index

            vertices = reshape([vertices, x], [2, size(vertices, 2) + 1])
            index = size(vertices, 2)
        end subroutine add_vertex

        subroutine add_edge(a, b, c, ua, ub)
            integer, intent(in) :: a, b, c
            real(dp), intent(in) :: ua, ub
            type(edge_t), allocatable :: grown(:)

            allocate(grown(size(edges) + 1))
            grown(:size(edges)) = edges
            grown(size(grown)) = edge_t(a, b, c, ua, ub)
            call move_alloc(grown, edges)
        end subroutine add_edge

        subroutine add_block(c)
            integer, intent(in) :: c(4)
            real(dp) :: area
            integer :: i0, i1

            area = 0.0_dp
            do i0 = 1, 4
                i1 = modulo(i0, 4) + 1
                area = area + vertices(1, c(i0))*vertices(2, c(i1)) - &
                    vertices(1, c(i1))*vertices(2, c(i0))
            end do
            nb = nb + 1
            if (area >= 0.0_dp) then
                corners(:, nb) = c
            else
                corners(:, nb) = [c(1), c(4), c(3), c(2)]
            end if
        end subroutine add_block

        subroutine side_window(a, b, c, u0, du)
            integer, intent(in) :: a, b
            integer, intent(out) :: c
            real(dp), intent(out) :: u0, du
            integer :: e

            c = 0
            u0 = 0.0_dp
            du = 0.0_dp
            do e = 1, size(edges)
                if (edges(e)%curve == 0) cycle
                if (edges(e)%a == a .and. edges(e)%b == b) then
                    c = edges(e)%curve
                    u0 = edges(e)%ua
                    du = edges(e)%ub - edges(e)%ua
                    return
                else if (edges(e)%a == b .and. edges(e)%b == a) then
                    c = edges(e)%curve
                    u0 = edges(e)%ub
                    du = edges(e)%ua - edges(e)%ub
                    return
                end if
            end do
        end subroutine side_window

        subroutine leg(samples, u, index)
            !! Leg curve from the X-point to its strike point on the wall.
            real(dp), intent(in) :: samples(:, :)
            real(dp), intent(out) :: u
            integer, intent(out) :: index
            real(dp), allocatable :: pts(:, :)
            real(dp) :: y(2), t(2)

            pts = samples
            pts(:, 1) = vertices(:, k(q))
            call closest_on_wall(curve(iw), pts(:, size(pts, 2)), u)
            call layout_curve_point(curve(iw), u, y, t)
            pts(:, size(pts, 2)) = y
            nc = nc + 1
            index = nc
            call open_curve(pts, p, LAYOUT_SEPARATRIX, curve(nc), status)
        end subroutine leg

    end subroutine layout_xpoint

    subroutine closed_curve(samples, p, tag, curve, status)
        !! Periodic spline of degree p through the samples at chord length.
        real(dp), intent(in) :: samples(:, :)
        integer, intent(in) :: p, tag
        type(layout_curve_t), intent(out) :: curve
        integer, intent(out) :: status

        real(dp), allocatable :: u(:), matrix(:, :), rhs(:, :), extraction(:, :, :)
        real(dp) :: value(8)
        integer, allocatable :: first(:), pivot(:)
        integer :: n, q, a, i

        status = 1
        n = size(samples, 2)
        if (n < p + 1) return
        allocate(u(n + 1))
        u(1) = 0.0_dp
        do q = 1, n
            u(q + 1) = u(q) + norm2(samples(:, modulo(q, n) + 1) - samples(:, q))
        end do
        u = u/u(n + 1)
        call build_periodic_extraction_1d(p, u, curve%knots, extraction, first, status)
        if (status /= 0) return
        allocate(matrix(n, n), rhs(n, 2), pivot(n))
        matrix = 0.0_dp
        do q = 1, n
            call local_bsplines(curve%knots, p, p + q, u(q), value(:p + 1))
            do a = 1, p + 1
                i = modulo(q + a - 2, n) + 1
                matrix(q, i) = matrix(q, i) + value(a)
            end do
            rhs(q, :) = samples(:, q)
        end do
        call dgesv(n, 2, matrix, n, pivot, rhs, n, status)
        if (status /= 0) return
        curve%tag = tag
        curve%degree = p
        curve%closed = .true.
        curve%offset = 0.0_dp
        curve%control = transpose(rhs)
    end subroutine closed_curve

    subroutine open_curve(samples, p, tag, curve, status)
        !! Clamped spline of degree min(p, n - 1) through the samples at chord
        !! length, knots by averaging (Piegl, Tiller, eq. 9.8).
        real(dp), intent(in) :: samples(:, :)
        integer, intent(in) :: p, tag
        type(layout_curve_t), intent(out) :: curve
        integer, intent(out) :: status

        real(dp), allocatable :: u(:), matrix(:, :), rhs(:, :)
        real(dp) :: value(8)
        integer, allocatable :: pivot(:)
        integer :: n, q, a, e, deg, nel, j

        status = 1
        n = size(samples, 2)
        if (size(samples, 1) /= 2 .or. n < 2) return
        deg = min(p, n - 1)
        allocate(u(n))
        u(1) = 0.0_dp
        do q = 2, n
            u(q) = u(q - 1) + norm2(samples(:, q) - samples(:, q - 1))
        end do
        if (u(n) <= 0.0_dp) return
        u = u/u(n)
        allocate(curve%knots(n + deg + 1))
        curve%knots(:deg + 1) = 0.0_dp
        curve%knots(n + 1:) = 1.0_dp
        do j = 1, n - deg - 1
            curve%knots(deg + 1 + j) = sum(u(j + 1:j + deg))/deg
        end do
        nel = n - deg
        allocate(matrix(n, n), rhs(n, 2), pivot(n))
        matrix = 0.0_dp
        do q = 1, n
            e = nel
            do a = 1, nel - 1
                if (u(q) < curve%knots(deg + a + 1)) then
                    e = a
                    exit
                end if
            end do
            call local_bsplines(curve%knots, deg, deg + e, u(q), value(:deg + 1))
            matrix(q, e:e + deg) = value(:deg + 1)
            rhs(q, :) = samples(:, q)
        end do
        call dgesv(n, 2, matrix, n, pivot, rhs, n, status)
        if (status /= 0) return
        curve%tag = tag
        curve%degree = deg
        curve%closed = .false.
        curve%offset = 0.0_dp
        curve%control = transpose(rhs)
    end subroutine open_curve

    subroutine closest_on_wall(wall, x, u)
        !! Parameter of the wall point closest to x (sampling, then Newton).
        type(layout_curve_t), intent(in) :: wall
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: u

        real(dp) :: y(2), t(2), best, f, h, y2(2), t2(2)
        integer :: q, it
        integer, parameter :: n = 2000

        best = huge(1.0_dp)
        u = 0.0_dp
        do q = 0, n - 1
            call layout_curve_point(wall, real(q, dp)/n, y, t)
            if (norm2(y - x) < best) then
                best = norm2(y - x)
                u = real(q, dp)/n
            end if
        end do
        h = 1.0e-6_dp
        do it = 1, 30
            call layout_curve_point(wall, u, y, t)
            f = dot_product(y - x, t)
            call layout_curve_point(wall, u + h, y2, t2)
            f = f/((dot_product(y2 - x, t2) - f)/h)
            u = modulo(u - f, 1.0_dp)
            if (abs(f) < 1.0e-15_dp) exit
        end do
    end subroutine closest_on_wall

    subroutine ray_to_wall(wall, center, x, u, status)
        !! Wall parameter where the ray from center through x meets the wall.
        type(layout_curve_t), intent(in) :: wall
        real(dp), intent(in) :: center(2), x(2)
        real(dp), intent(out) :: u
        integer, intent(out) :: status

        real(dp) :: dir(2), y(2), t(2), f0, f1, a, b, c, fc
        integer :: q, it
        integer, parameter :: n = 2000

        status = 2
        dir = x - center
        u = 0.0_dp
        do q = 0, n - 1
            a = real(q, dp)/n
            b = real(q + 1, dp)/n
            call layout_curve_point(wall, a, y, t)
            f0 = cross(y - center, dir)
            if (dot_product(y - center, dir) <= 0.0_dp) cycle
            call layout_curve_point(wall, b, y, t)
            f1 = cross(y - center, dir)
            if (f0*f1 > 0.0_dp) cycle
            do it = 1, 60
                c = 0.5_dp*(a + b)
                call layout_curve_point(wall, c, y, t)
                fc = cross(y - center, dir)
                if (f0*fc <= 0.0_dp) then
                    b = c
                else
                    a = c
                    f0 = fc
                end if
            end do
            u = 0.5_dp*(a + b)
            status = 0
            return
        end do

    contains

        pure real(dp) function cross(v, w)
            real(dp), intent(in) :: v(2), w(2)
            cross = v(1)*w(2) - v(2)*w(1)
        end function cross

    end subroutine ray_to_wall

end module fortfem_spline_xpoint
