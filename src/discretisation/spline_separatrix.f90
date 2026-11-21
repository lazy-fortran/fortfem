module fortfem_spline_separatrix
    !! Block layouts aligned with a separatrix through an X-point inside the
    !! domain (doc/design/multipatch_spline_spaces.md, separatrix layouts).
    !!
    !! Doublet: two lobes meet at the X-point X. Each lobe is bounded by a
    !! loop from X counter-clockwise about its O-point back to X. The lobe is
    !! a 5-block O-grid whose outer corners are X and the loop points at
    !! chord-length parameter 1/4, 1/2, 3/4, and whose radial lines run
    !! along rays from the O-point. The region between the two loops and the
    !! wall is a ring of 8 blocks: their inner sides are the 8 loop pieces,
    !! the radial lines continue the lobe rays to the wall, and at X two
    !! lines run to the wall along the bisectors of the two outer sectors.
    !! X is a vertex of valence 8. Every loop piece carries LAYOUT_SEPARATRIX
    !! on both of its blocks, and every block map is the Coons patch of its
    !! four sides.
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

    type :: edge_t
        integer :: a = 0, b = 0, curve = 0
        real(dp) :: ua = 0.0_dp, ub = 0.0_dp
    end type edge_t

    public :: layout_doublet

contains

    subroutine layout_doublet(curves, core_fraction, layout, status, degree)
        !! Doublet layout of 18 blocks. The topology does not depend on the
        !! samples, so moved curves keep the spline extraction
        !! (spline_space_update_geometry). Status 1: invalid input, 2: the
        !! lobes or the wall rays are not in counter-clockwise order about X,
        !! or a ray misses the wall, 3: a block is not counter-clockwise or
        !! the layout is not conforming.
        type(doublet_curves_t), intent(in) :: curves
        real(dp), intent(in) :: core_fraction
        type(block_layout_t), intent(out) :: layout
        integer, intent(out) :: status
        integer, intent(in), optional :: degree

        real(dp), parameter :: pi = acos(-1.0_dp)
        type(layout_curve_t) :: curve(3)
        type(edge_t), allocatable :: edges(:)
        real(dp), allocatable :: vertices(:, :), points(:, :)
        integer, allocatable :: corners(:, :), tags(:, :)
        real(dp) :: x(2), tangent(2), start(2), finish(2), uwall(8), turn
        real(dp) :: bisector(2, 2), origin(2), dir(2), angle_start(2), angle_end(2)
        integer :: p, k, j, jn, nb, ix, outer(5, 2), inner(4, 2), ring(8), w(8)
        integer :: s, i

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
        call closed_curve(curves%wall, p, LAYOUT_WALL, curve(1), status)
        if (status /= 0) return
        do k = 1, 2
            points = curves%lobe(k)%points
            points(:, 1) = curves%xpoint
            points(:, size(points, 2)) = curves%xpoint
            call open_curve(points, p, LAYOUT_SEPARATRIX, curve(1 + k), status)
            if (status /= 0) return
        end do
        status = 2
        do k = 1, 2
            if (loop_area(curves%lobe(k)%points) <= 0.0_dp) return
            call layout_curve_point(curve(1 + k), 0.0_dp, x, start)
            call layout_curve_point(curve(1 + k), 1.0_dp, x, finish)
            angle_start(k) = atan2(start(2), start(1))
            angle_end(k) = atan2(-finish(2), -finish(1))
        end do
        turn = sector(angle_start(1), angle_end(1)) + sector(angle_end(1), &
            angle_start(2)) + sector(angle_start(2), angle_end(2)) + &
            sector(angle_end(2), angle_start(1))
        if (abs(turn - 2.0_dp*pi) > 1.0e-8_dp) return
        bisector(:, 1) = direction(angle_end(2) + 0.5_dp*sector(angle_end(2), &
            angle_start(1)))
        bisector(:, 2) = direction(angle_end(1) + 0.5_dp*sector(angle_end(1), &
            angle_start(2)))

        allocate(vertices(2, 0), edges(0))
        call add_vertex(curves%xpoint, ix)
        do k = 1, 2
            outer(1, k) = ix
            outer(5, k) = ix
            do j = 2, 4
                call layout_curve_point(curve(1 + k), 0.25_dp*(j - 1), x, tangent)
                call add_vertex(x, outer(j, k))
            end do
            do j = 1, 4
                call add_vertex(curves%opoint(:, k) + core_fraction* &
                    (vertices(:, outer(j, k)) - curves%opoint(:, k)), inner(j, k))
            end do
            do j = 1, 4
                call add_edge(outer(j, k), outer(j + 1, k), 1 + k, 0.25_dp*(j - 1), &
                    0.25_dp*j)
            end do
        end do
        ! the ring runs counter-clockwise about X: the outer sector before
        ! lobe 1, lobe 1, the outer sector before lobe 2, lobe 2
        ring = [ix, outer(2:4, 1), ix, outer(2:4, 2)]
        do j = 1, 8
            k = merge(1, 2, j <= 4)
            if (j == 1 .or. j == 5) then
                origin = curves%xpoint
                dir = bisector(:, merge(1, 2, j == 1))
            else
                origin = curves%opoint(:, k)
                dir = vertices(:, ring(j)) - origin
            end if
            call ray_to_wall(curve(1), origin, origin + dir, uwall(j), status)
            if (status /= 0) return
            call layout_curve_point(curve(1), uwall(j), x, tangent)
            call add_vertex(x, w(j))
        end do
        status = 2
        turn = 0.0_dp
        do j = 1, 8
            jn = modulo(j, 8) + 1
            if (modulo(uwall(jn) - uwall(j), 1.0_dp) <= 0.0_dp) return
            turn = turn + modulo(uwall(jn) - uwall(j), 1.0_dp)
        end do
        if (abs(turn - 1.0_dp) > 1.0e-8_dp) return
        do j = 1, 8
            jn = modulo(j, 8) + 1
            call add_edge(w(j), w(jn), 1, uwall(j), uwall(j) + &
                modulo(uwall(jn) - uwall(j), 1.0_dp))
        end do

        allocate(corners(4, 18))
        nb = 0
        do k = 1, 2
            call add_block(inner(:, k))
            do j = 1, 4
                jn = modulo(j, 4) + 1
                call add_block([outer(j, k), outer(j + 1, k), inner(jn, k), &
                    inner(j, k)])
            end do
        end do
        do j = 1, 8
            jn = modulo(j, 8) + 1
            call add_block([ring(j), w(j), w(jn), ring(jn)])
        end do
        allocate(tags(4, nb))
        tags = LAYOUT_WALL
        status = 3
        call layout_from_blocks(vertices, corners(:, :nb), tags, layout, i)
        if (i /= 0) return
        layout%curves = curve
        layout%coons = .true.
        layout%center = curves%xpoint
        do j = 1, nb
            do s = 1, 4
                call side_window(corners(side_start(s), j), corners(side_end(s), j), &
                    layout%side_curve(s, j), layout%side_u0(s, j), &
                    layout%side_du(s, j))
                if (layout%side_curve(s, j) >= 2) &
                    layout%side_tag(s, j) = LAYOUT_SEPARATRIX
            end do
        end do
        status = 0

    contains

        subroutine add_vertex(y, index)
            real(dp), intent(in) :: y(2)
            integer, intent(out) :: index

            vertices = reshape([vertices, y], [2, size(vertices, 2) + 1])
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
            !! Corners in either orientation, stored counter-clockwise.
            integer, intent(in) :: c(4)

            nb = nb + 1
            if (loop_area(vertices(:, c)) >= 0.0_dp) then
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

    end subroutine layout_doublet

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
