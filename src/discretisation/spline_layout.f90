module fortfem_spline_layout
    !! Quadrilateral block layouts of a poloidal reference plane for multipatch
    !! spline spaces (doc/design/multipatch_spline_spaces.md).
    !!
    !! Block corners are listed at xi = (0,0), (1,0), (1,1), (0,1). Sides are
    !! 1: xi2 = 0, 2: xi1 = 1, 3: xi2 = 1, 4: xi1 = 0, each oriented along its
    !! increasing parameter. Interior sides are straight. A boundary side may be
    !! a circular arc, and the block map is then ruled between that arc and the
    !! opposite straight side.
    use fortfem_kinds, only: dp
    implicit none
    private

    integer, parameter, public :: LAYOUT_INTERIOR = 0
    integer, parameter, public :: LAYOUT_BOUNDARY = 1
    integer, parameter, public :: LAYOUT_WALL = 2
    integer, parameter, public :: LAYOUT_PLASMA_BOUNDARY = 3
    integer, parameter, public :: LAYOUT_SEPARATRIX = 4

    integer, parameter, public :: side_start(4) = [1, 2, 4, 1]
    integer, parameter, public :: side_end(4) = [2, 3, 3, 4]
    integer, parameter, public :: opposite_side(4) = [3, 4, 1, 2]

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
    end type block_layout_t

    public :: layout_from_blocks
    public :: layout_set_arc
    public :: layout_ogrid
    public :: layout_side_point

contains

    subroutine layout_from_blocks(vertices, corners, boundary_tag, layout, status)
        !! Connect blocks that share a side (same two vertices). Status 2: a
        !! vertex pair belongs to more than two sides, or a block is not
        !! counter-clockwise.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: corners(:, :), boundary_tag(:, :)
        type(block_layout_t), intent(out) :: layout
        integer, intent(out) :: status

        integer :: nb, nv, k, s, k2, s2, a, b, a2, b2

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
            layout%arc_side(nb), layout%arc_center(2, nb))
        layout%neighbour = 0
        layout%neighbour_side = 0
        layout%reversed = .false.
        layout%arc_side = 0
        layout%arc_center = 0.0_dp
        status = 2
        do k = 1, nb
            if (signed_area(vertices(:, corners(:, k))) <= 0.0_dp) return
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

    subroutine layout_ogrid(center, radius, core_fraction, layout, status)
        !! Five-block O-grid of the disc |s - center| <= radius: a centre
        !! square of half-width core_fraction*radius/sqrt(2) and four ruled
        !! blocks with arc outer sides. The four inner corners are valence-3
        !! vertices at mid-radius. Outer sides carry LAYOUT_WALL.
        real(dp), intent(in) :: center(2), radius, core_fraction
        type(block_layout_t), intent(out) :: layout
        integer, intent(out) :: status

        real(dp) :: vertices(2, 8), a, r
        integer :: corners(4, 5), tags(4, 5), k

        status = 1
        if (radius <= 0.0_dp .or. core_fraction <= 0.0_dp .or. &
            core_fraction >= 1.0_dp) return
        a = core_fraction*radius/sqrt(2.0_dp)
        r = radius/sqrt(2.0_dp)
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
        tags = LAYOUT_WALL
        call layout_from_blocks(vertices, corners, tags, layout, status)
        if (status /= 0) return
        do k = 2, 5
            call layout_set_arc(layout, k, 1, center, status)
            if (status /= 0) return
        end do
    end subroutine layout_ogrid

    pure subroutine layout_side_point(layout, block, side, t, point, tangent)
        !! Point and d/dt of a block side at parameter t in [0, 1].
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: block, side
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2), tangent(2)

        real(dp) :: a(2), b(2), c(2), angle_a, angle_b, angle, radius, pi

        a = layout%vertices(:, layout%corners(side_start(side), block))
        b = layout%vertices(:, layout%corners(side_end(side), block))
        if (layout%arc_side(block) /= side) then
            point = a + t*(b - a)
            tangent = b - a
            return
        end if
        pi = acos(-1.0_dp)
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
