module fortfem_spline_layout
    !! Quadrilateral block layouts of a poloidal reference plane for multipatch
    !! spline spaces (doc/design/multipatch_spline_spaces.md).
    !!
    !! Block corners are listed at xi = (0,0), (1,0), (1,1), (0,1). Sides are
    !! 1: xi2 = 0, 2: xi1 = 1, 3: xi2 = 1, 4: xi1 = 0, each oriented along its
    !! increasing parameter. Interior sides are straight. A boundary side may be
    !! a circular arc or a piece of a closed curve, and the block map is then
    !! ruled between that side and the opposite straight side.
    !!
    !! Closed curves are periodic splines of odd degree through counter-
    !! clockwise samples at their normalised chord-length parameter u, u = 0
    !! at the first sample. A curve side runs over u = offset + u0 + du t,
    !! with u0 and du fixed by the topology (layout_ogrid: the core corners at
    !! u = offset + 0, 1/4, 1/2, 3/4, offset = 5/8, so the first sample on the
    !! outboard midplane puts the corners near the diagonals).
    use fortfem_kinds, only: dp
    use fortfem_spline_extraction_1d, only: build_periodic_extraction_1d, &
        local_bsplines
    implicit none
    private

    integer, parameter, public :: LAYOUT_INTERIOR = 0
    integer, parameter, public :: LAYOUT_BOUNDARY = 1
    integer, parameter, public :: LAYOUT_WALL = 2
    integer, parameter, public :: LAYOUT_PLASMA_BOUNDARY = 3
    integer, parameter, public :: LAYOUT_SEPARATRIX = 4
    !> Sides where a periodic annulus meets the core blocks (C0, the core
    !> traces lie in the annulus trace space).
    integer, parameter, public :: LAYOUT_ANNULUS_INTERFACE = 5

    integer, parameter, public :: side_start(4) = [1, 2, 4, 1]
    integer, parameter, public :: side_end(4) = [2, 3, 3, 4]
    integer, parameter, public :: opposite_side(4) = [3, 4, 1, 2]

    type, public :: layout_curve_t
        integer :: tag = 0
        integer :: degree = 3
        logical :: closed = .true. ! else clamped on [0, 1]
        real(dp) :: offset = 0.625_dp
        real(dp), allocatable :: knots(:) ! periodically extended
        real(dp), allocatable :: control(:, :) ! (2, n) periodic control points
    end type layout_curve_t

    type, public :: block_layout_t
        integer :: block_count = 0
        integer :: vertex_count = 0
        real(dp), allocatable :: vertices(:, :) ! (2, vertex)
        integer, allocatable :: corners(:, :) ! (4, block)
        integer, allocatable :: side_tag(:, :) ! (4, block)
        integer, allocatable :: neighbour(:, :) ! (4, block), 0 on the boundary
        integer, allocatable :: neighbour_side(:, :) ! (4, block)
        logical, allocatable :: reversed(:, :) ! (4, block)
        integer, allocatable :: valence(:) ! blocks at each vertex
        logical, allocatable :: boundary_vertex(:)
        integer, allocatable :: arc_side(:) ! (block), 0 for none
        real(dp), allocatable :: arc_center(:, :) ! (2, block)
        !> Periodic in xi1 (sides 2 and 4 coincide): an annulus between the
        !> full circles of sides 1 and 3 about arc_center.
        logical, allocatable :: periodic(:) ! (block)
        type(layout_curve_t), allocatable :: curves(:)
        !> Interior sides coupled only in C0 in a C1 space (both sides set).
        logical, allocatable :: side_c0(:, :) ! (4, block)
        !> Block map is the Coons patch of its four sides (any may be curved).
        logical, allocatable :: coons(:) ! (block)
        !> Block directions refined by radial (anisotropic) THB levels.
        logical, allocatable :: radial_refine(:, :) ! (2, block)
        integer, allocatable :: side_curve(:, :) ! (4, block), 0 for none
        real(dp), allocatable :: side_u0(:, :), side_du(:, :) ! (4, block)
        !> O-grid core: inner corner k = center + core_scale (outer corner
        !> k + 4 - center), kept when the outer curve changes.
        real(dp) :: core_scale = 0.0_dp
        real(dp) :: center(2) = 0.0_dp
    end type block_layout_t

    public :: layout_from_blocks
    public :: layout_set_arc
    public :: layout_ogrid
    public :: layout_side_point
    public :: layout_set_curve
    public :: layout_curve_point
    public :: layout_set_c0

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

contains

    subroutine layout_from_blocks(vertices, corners, boundary_tag, layout, status, &
            oriented)
        !! Connect blocks that share a side (same two vertices). Status 2: a
        !! vertex pair belongs to more than two sides, or a block is not
        !! counter-clockwise. oriented = .true.: the caller has checked the
        !! orientation along curved sides, where a thin curved block can have
        !! a clockwise corner polygon.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: corners(:, :), boundary_tag(:, :)
        type(block_layout_t), intent(out) :: layout
        integer, intent(out) :: status
        logical, intent(in), optional :: oriented

        integer :: nb, nv, k, s, k2, s2, a, b, a2, b2
        logical :: check

        status = 1
        nb = size(corners, 2)
        nv = size(vertices, 2)
        if (size(vertices, 1) /= 2 .or. size(corners, 1) /= 4 .or. nb < 1) return
        if (any(shape(boundary_tag) /= [4, nb])) return
        if (any(corners < 1) .or. any(corners > nv)) return
        layout%block_count = nb
        layout%vertex_count = nv
        layout%vertices = vertices
        layout%corners = corners
        allocate(layout%side_tag(4, nb), layout%neighbour(4, nb), &
            layout%neighbour_side(4, nb), layout%reversed(4, nb), &
            layout%valence(nv), layout%boundary_vertex(nv), &
            layout%arc_side(nb), layout%arc_center(2, nb), layout%periodic(nb), &
            layout%side_curve(4, nb), layout%side_u0(4, nb), layout%side_du(4, nb), &
            layout%curves(0), layout%side_c0(4, nb), layout%coons(nb), &
            layout%radial_refine(2, nb))
        layout%radial_refine = .false.
        layout%side_c0 = .false.
        layout%coons = .false.
        layout%periodic = .false.
        layout%side_curve = 0
        layout%side_u0 = 0.0_dp
        layout%side_du = 0.0_dp
        layout%neighbour = 0
        layout%neighbour_side = 0
        layout%reversed = .false.
        layout%arc_side = 0
        layout%arc_center = 0.0_dp
        status = 2
        check = .true.
        if (present(oriented)) check = .not. oriented
        do k = 1, nb
            if (check .and. signed_area(vertices(:, corners(:, k))) <= 0.0_dp) return
        end do
        do k = 1, nb
            do s = 1, 4
                a = corners(side_start(s), k)
                b = corners(side_end(s), k)
                do k2 = 1, nb
                    do s2 = 1, 4
                        if (k2 == k .and. s2 == s) cycle
                        a2 = corners(side_start(s2), k2)
                        b2 = corners(side_end(s2), k2)
                        if (.not. ((a == a2 .and. b == b2) .or. &
                            (a == b2 .and. b == a2))) cycle
                        if (layout%neighbour(s, k) /= 0) return
                        layout%neighbour(s, k) = k2
                        layout%neighbour_side(s, k) = s2
                        layout%reversed(s, k) = a /= a2
                    end do
                end do
            end do
        end do
        where (layout%neighbour == 0)
            layout%side_tag = max(boundary_tag, LAYOUT_BOUNDARY)
        elsewhere
            layout%side_tag = LAYOUT_INTERIOR
        end where
        layout%valence = 0
        layout%boundary_vertex = .false.
        do k = 1, nb
            do s = 1, 4
                layout%valence(corners(s, k)) = layout%valence(corners(s, k)) + 1
                if (layout%neighbour(s, k) == 0) then
                    layout%boundary_vertex(corners(side_start(s), k)) = .true.
                    layout%boundary_vertex(corners(side_end(s), k)) = .true.
                end if
            end do
        end do
        status = 0
    end subroutine layout_from_blocks

    subroutine layout_set_arc(layout, block, side, center, status)
        !! Make a boundary side the circular arc about center through its two
        !! end vertices (both at the same distance from center).
        type(block_layout_t), intent(inout) :: layout
        integer, intent(in) :: block, side
        real(dp), intent(in) :: center(2)
        integer, intent(out) :: status

        real(dp) :: ra, rb

        status = 1
        if (block < 1 .or. block > layout%block_count .or. side < 1 .or. side > 4) &
            return
        if (layout%neighbour(side, block) /= 0) return
        ra = norm2(layout%vertices(:, layout%corners(side_start(side), block)) - &
            center)
        rb = norm2(layout%vertices(:, layout%corners(side_end(side), block)) - &
            center)
        if (abs(ra - rb) > 1.0e-12_dp*max(ra, rb) .or. ra <= 0.0_dp) return
        layout%arc_side(block) = side
        layout%arc_center(:, block) = center
        status = 0
    end subroutine layout_set_arc

    subroutine layout_ogrid(center, radius, core_fraction, layout, status, polygon, &
            interface_fraction)
        !! Five-block O-grid of the disc |s - center| <= radius: a centre
        !! square of half-width core_fraction*radius/sqrt(2) and four ruled
        !! blocks with arc outer sides. The four inner corners are valence-3
        !! vertices at mid-radius. Outer sides carry LAYOUT_WALL. With
        !! polygon = .true. the outer sides stay straight (the square
        !! inscribed in the disc), so every block is bilinear.
        !! With interface_fraction the five blocks fill the disc of radius
        !! interface_fraction*radius, and block 6 is the periodic annulus up
        !! to the wall: xi1 runs counter-clockwise from the angle of vertex 5,
        !! xi2 = 0 is the wall (side 1) and xi2 = 1 the interface (side 3).
        real(dp), intent(in) :: center(2), radius, core_fraction
        type(block_layout_t), intent(out) :: layout
        integer, intent(out) :: status
        logical, intent(in), optional :: polygon
        real(dp), intent(in), optional :: interface_fraction

        real(dp) :: vertices(2, 8), a, r, outer
        integer :: corners(4, 5), tags(4, 5), k
        logical :: straight

        status = 1
        straight = .false.
        if (present(polygon)) straight = polygon
        outer = 1.0_dp
        if (present(interface_fraction)) then
            if (straight .or. interface_fraction <= core_fraction .or. &
                interface_fraction >= 1.0_dp) return
            outer = interface_fraction
        end if
        if (radius <= 0.0_dp .or. core_fraction <= 0.0_dp .or. &
            core_fraction >= 1.0_dp) return
        a = core_fraction*radius/sqrt(2.0_dp)
        r = outer*radius/sqrt(2.0_dp)
        vertices(:, 1) = center + [-a, -a]
        vertices(:, 2) = center + [a, -a]
        vertices(:, 3) = center + [a, a]
        vertices(:, 4) = center + [-a, a]
        vertices(:, 5) = center + [-r, -r]
        vertices(:, 6) = center + [r, -r]
        vertices(:, 7) = center + [r, r]
        vertices(:, 8) = center + [-r, r]
        corners(:, 1) = [1, 2, 3, 4]
        corners(:, 2) = [5, 6, 2, 1]
        corners(:, 3) = [6, 7, 3, 2]
        corners(:, 4) = [7, 8, 4, 3]
        corners(:, 5) = [8, 5, 1, 4]
        tags = merge(LAYOUT_ANNULUS_INTERFACE, LAYOUT_WALL, &
            present(interface_fraction))
        call layout_from_blocks(vertices, corners, tags, layout, status)
        if (status /= 0) return
        layout%core_scale = core_fraction/outer
        layout%center = center
        do k = 2, 5
            layout%side_u0(1, k) = 0.25_dp*(k - 2)
            layout%side_du(1, k) = 0.25_dp
        end do
        if (straight) return
        ! C1 cannot reach a curved side at a core corner (doc, Model)
        layout%side_c0(2:4:2, 2:5) = .true.
        do k = 2, 5
            call layout_set_arc(layout, k, 1, center, status)
            if (status /= 0) return
        end do
        if (present(interface_fraction)) call append_annulus(layout, center, &
            center + (radius/outer)*(vertices(:, 5) - center), 5)
    end subroutine layout_ogrid

    subroutine append_annulus(layout, center, wall_point, interface_vertex)
        !! Periodic block 6 with corners (w, w, v, v): w on the wall at the
        !! angle of the interface vertex v.
        type(block_layout_t), intent(inout) :: layout
        real(dp), intent(in) :: center(2), wall_point(2)
        integer, intent(in) :: interface_vertex

        integer :: nb, nv

        nb = layout%block_count + 1
        nv = layout%vertex_count + 1
        layout%block_count = nb
        layout%vertex_count = nv
        layout%vertices = reshape([layout%vertices, wall_point], [2, nv])
        layout%corners = reshape([layout%corners, nv, nv, interface_vertex, &
            interface_vertex], [4, nb])
        layout%side_tag = reshape([layout%side_tag, LAYOUT_WALL, LAYOUT_INTERIOR, &
            LAYOUT_ANNULUS_INTERFACE, LAYOUT_INTERIOR], [4, nb])
        layout%neighbour = reshape([layout%neighbour, 0, 0, 0, 0], [4, nb])
        layout%neighbour_side = reshape([layout%neighbour_side, 0, 0, 0, 0], [4, nb])
        layout%reversed = reshape([layout%reversed, .false., .false., .false., &
            .false.], [4, nb])
        layout%valence = [layout%valence, 1]
        layout%valence(interface_vertex) = layout%valence(interface_vertex) + 1
        layout%boundary_vertex = [layout%boundary_vertex, .true.]
        layout%arc_side = [layout%arc_side, 1]
        layout%arc_center = reshape([layout%arc_center, center], [2, nb])
        layout%periodic = [layout%periodic, .true.]
        layout%side_curve = reshape([layout%side_curve, 0, 0, 0, 0], [4, nb])
        layout%side_c0 = reshape([layout%side_c0, .false., .false., .false., &
            .false.], [4, nb])
        layout%coons = [layout%coons, .false.]
        layout%radial_refine = reshape([layout%radial_refine, .false., .true.], &
            [2, nb])
        layout%side_u0 = reshape([layout%side_u0, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
            [4, nb])
        layout%side_du = reshape([layout%side_du, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], &
            [4, nb])
    end subroutine append_annulus

    subroutine layout_set_c0(layout, block, side, c0, status)
        !! Couple an interior side (and its neighbour side) only in C0 in C1
        !! spaces. Status 1: invalid block or side, or a boundary side.
        type(block_layout_t), intent(inout) :: layout
        integer, intent(in) :: block, side
        logical, intent(in) :: c0
        integer, intent(out) :: status

        integer :: k2

        status = 1
        if (block < 1 .or. block > layout%block_count .or. side < 1 .or. side > 4) &
            return
        k2 = layout%neighbour(side, block)
        if (k2 == 0) return
        layout%side_c0(side, block) = c0
        layout%side_c0(layout%neighbour_side(side, block), k2) = c0
        status = 0
    end subroutine layout_set_c0

    subroutine layout_set_curve(layout, tag, samples, status, degree)
        !! Make every side with the given tag a piece of the closed curve
        !! through samples(2, n) (counter-clockwise, first sample on the
        !! outboard midplane), and move the side end vertices onto it. On an
        !! O-grid whose core sides carry the tag, the inner corners follow.
        !! The topology, and with it the spline extraction, is unchanged.
        !! Status 1: invalid input or no curve side with this tag,
        !! 2: the samples are clockwise or a block turns over.
        type(block_layout_t), intent(inout) :: layout
        integer, intent(in) :: tag
        real(dp), intent(in) :: samples(:, :)
        integer, intent(out) :: status
        integer, intent(in), optional :: degree

        type(layout_curve_t) :: curve
        type(layout_curve_t), allocatable :: curves(:)
        real(dp), allocatable :: u(:), matrix(:, :), rhs(:, :), extraction(:, :, :)
        real(dp) :: value(8), tangent(2), length, area
        integer, allocatable :: first(:), pivot(:)
        integer :: n, p, q, i, k, s, index, a

        status = 1
        p = 3
        if (present(degree)) p = degree
        n = size(samples, 2)
        if (size(samples, 1) /= 2 .or. p < 1 .or. p > 7 .or. mod(p, 2) /= 1 .or. &
            n < p + 1) return
        if (.not. any(layout%side_tag == tag .and. layout%side_du /= 0.0_dp)) return
        area = 0.0_dp
        do q = 1, n
            i = modulo(q, n) + 1
            area = area + samples(1, q)*samples(2, i) - samples(1, i)*samples(2, q)
        end do
        status = 2
        if (area <= 0.0_dp) return
        allocate(u(n + 1))
        u(1) = 0.0_dp
        do q = 1, n
            u(q + 1) = u(q) + norm2(samples(:, modulo(q, n) + 1) - samples(:, q))
        end do
        length = u(n + 1)
        status = 1
        if (length <= 0.0_dp) return
        u = u/length
        u(n + 1) = 1.0_dp
        call build_periodic_extraction_1d(p, u, curve%knots, extraction, first, &
            status)
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
        if (status /= 0) then
            status = 1
            return
        end if
        curve%tag = tag
        curve%degree = p
        curve%control = transpose(rhs)
        index = 0
        do i = 1, size(layout%curves)
            if (layout%curves(i)%tag == tag) index = i
        end do
        if (index == 0) then
            allocate(curves(size(layout%curves) + 1))
            curves(:size(layout%curves)) = layout%curves
            index = size(curves)
            call move_alloc(curves, layout%curves)
        end if
        layout%curves(index) = curve
        do k = 1, layout%block_count
            do s = 1, 4
                if (layout%side_tag(s, k) /= tag .or. layout%side_du(s, k) == 0.0_dp) &
                    cycle
                layout%side_curve(s, k) = index
                if (.not. layout%periodic(k)) layout%arc_side(k) = s
                call layout_curve_point(curve, curve%offset + layout%side_u0(s, k), &
                    layout%vertices(:, layout%corners(side_start(s), k)), tangent)
                call layout_curve_point(curve, curve%offset + layout%side_u0(s, k) + &
                    layout%side_du(s, k), &
                    layout%vertices(:, layout%corners(side_end(s), k)), tangent)
            end do
        end do
        if (layout%core_scale > 0.0_dp .and. layout%block_count >= 5) then
            if (layout%side_tag(1, 2) == tag) then
                do k = 1, 4
                    layout%vertices(:, k) = layout%center + layout%core_scale* &
                        (layout%vertices(:, k + 4) - layout%center)
                end do
            end if
        end if
        status = 2
        do k = 1, layout%block_count
            if (layout%periodic(k)) cycle
            if (signed_area(layout%vertices(:, layout%corners(:, k))) <= 0.0_dp) &
                return
        end do
        status = 0
    end subroutine layout_set_curve

    pure subroutine layout_curve_point(curve, u, point, tangent)
        !! Point and d/du of a curve at u (closed: modulo 1, open: clamped).
        type(layout_curve_t), intent(in) :: curve
        real(dp), intent(in) :: u
        real(dp), intent(out) :: point(2), tangent(2)

        real(dp) :: x, value(curve%degree + 1), lower(curve%degree), d
        integer :: p, n, e, a, i, nel

        p = curve%degree
        n = size(curve%control, 2)
        nel = size(curve%knots) - 2*p - 1
        if (curve%closed) then
            x = modulo(u, 1.0_dp)
        else
            x = min(max(u, 0.0_dp), 1.0_dp)
        end if
        e = nel
        do a = 1, nel - 1
            if (x < curve%knots(p + a + 1)) then
                e = a
                exit
            end if
        end do
        call local_bsplines(curve%knots, p, p + e, x, value)
        call local_bsplines(curve%knots, p - 1, p + e, x, lower)
        point = 0.0_dp
        tangent = 0.0_dp
        do a = 1, p + 1
            i = e + a - 1
            point = point + value(a)*curve%control(:, modulo(i - 1, n) + 1)
            d = 0.0_dp
            if (a >= 2) d = d + lower(a - 1)/(curve%knots(i + p) - curve%knots(i))
            if (a <= p) d = d - lower(a)/(curve%knots(i + p + 1) - curve%knots(i + 1))
            tangent = tangent + p*d*curve%control(:, modulo(i - 1, n) + 1)
        end do
    end subroutine layout_curve_point

    pure subroutine layout_side_point(layout, block, side, t, point, tangent)
        !! Point and d/dt of a block side at parameter t in [0, 1].
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: block, side
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2), tangent(2)

        real(dp) :: a(2), b(2), c(2), angle_a, angle_b, angle, radius, pi

        if (layout%side_curve(side, block) /= 0) then
            associate (curve => layout%curves(layout%side_curve(side, block)))
                call layout_curve_point(curve, curve%offset + &
                    layout%side_u0(side, block) + layout%side_du(side, block)*t, &
                    point, tangent)
                tangent = layout%side_du(side, block)*tangent
            end associate
            return
        end if
        a = layout%vertices(:, layout%corners(side_start(side), block))
        b = layout%vertices(:, layout%corners(side_end(side), block))
        pi = acos(-1.0_dp)
        if (layout%periodic(block) .and. (side == 1 .or. side == 3)) then
            c = layout%arc_center(:, block)
            radius = norm2(a - c)
            angle = atan2(a(2) - c(2), a(1) - c(1)) + 2.0_dp*pi*t
            point = c + radius*[cos(angle), sin(angle)]
            tangent = 2.0_dp*pi*radius*[-sin(angle), cos(angle)]
            return
        end if
        if (layout%arc_side(block) /= side) then
            point = a + t*(b - a)
            tangent = b - a
            return
        end if
        c = layout%arc_center(:, block)
        radius = norm2(a - c)
        angle_a = atan2(a(2) - c(2), a(1) - c(1))
        angle_b = atan2(b(2) - c(2), b(1) - c(1))
        if (angle_b - angle_a > pi) angle_b = angle_b - 2.0_dp*pi
        if (angle_b - angle_a < -pi) angle_b = angle_b + 2.0_dp*pi
        angle = angle_a + t*(angle_b - angle_a)
        point = c + radius*[cos(angle), sin(angle)]
        tangent = radius*(angle_b - angle_a)*[-sin(angle), cos(angle)]
    end subroutine layout_side_point

    pure real(dp) function signed_area(p) result(area)
        real(dp), intent(in) :: p(2, 4)
        integer :: i, j

        area = 0.0_dp
        do i = 1, 4
            j = modulo(i, 4) + 1
            area = area + p(1, i)*p(2, j) - p(1, j)*p(2, i)
        end do
        area = 0.5_dp*area
    end function signed_area

end module fortfem_spline_layout
